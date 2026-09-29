# devtools::install_github("https://github.com/ErikKusch/KrigR", ref = "Development")

library(terra)
library(lubridate)

#------------------------------------------------------------------------------------------------------------
# Extract climate data from rasters and calculate rabbit climate variables
#------------------------------------------------------------------------------------------------------------


format_rabbit_climate <- function(data_dir = file.path("environmental_data", "data", "CDS"), 
                                  out_dir  = file.path("rabbit", "data", "model_input", "climate_historic"), 
                                  year_min = 2006, 
                                  year_max = 2025,
                                  template_file = file.path("environmental_data", "data", "template_raster_500.tif")){
  
  Dir.Data <- data_dir
  Dir.Out <- out_dir
  
  cat('directory for climate data: ', data_dir, fill = T)
  cat('directory for output files: ', out_dir, fill = T)
  
  
  ## create directories, if they don't exist yet
  if (!dir.exists(Dir.Out)) {dir.create(Dir.Out)}
  
  template <- rast(template_file)
  
  tas_files <- list.files(Dir.Data, pattern = "Kriged.nc", full.names = T)
  pr_file <- list.files(Dir.Data, pattern = "precipitation.grib", full.names = T)
  
  cat('Extracting Temperature')
  
  tas_rast <- rast(tas_files)
  
  tas_rast <- project(tas_rast, crs(template))
  tas_rast <- terra::resample(tas_rast, template)
  
  temp_donana <- as.data.frame(tas_rast, xy = T) 
  temp_donana[,-c(1:2)] <- temp_donana[,-c(1:2)] - 273.15
  
  temp_donana$x <- colFromX(tas_rast, temp_donana$x)
  temp_donana$y <- rowFromY(tas_rast, temp_donana$y)
  
  colnames(temp_donana) <- gsub("_Krigged_Kriged", "", colnames(temp_donana))
  y <- regmatches(colnames(temp_donana), regexpr("\\d{4}", colnames(temp_donana)))
  m <- regmatches(colnames(temp_donana), regexpr("\\d{1,2}$", colnames(temp_donana)))
  colnames(temp_donana)[-c(1:2)] <- paste(y, sprintf("%02d", as.integer(m)), sep = "_")
  
  cat('Extracting Precipitation')
  if(!file.exists(file.path("environmental_data", "data", "CDS", "precipitation.grib"))){ cat('precipitation data not found') }
  
  pr_rast <- rast(file.path(Dir.Data, "precipitation.grib"))
  pr_rast <- project(pr_rast, crs(template))
  pr_rast <- terra::resample(pr_rast, template)
  
  pr_donana <- as.data.frame(pr_rast, xy = T)
  pr_donana[,-c(1:2)] <- pr_donana[,-c(1:2)] * 1000 * 30.5
  
  pr_donana$x <- colFromX(pr_rast, pr_donana$x)
  pr_donana$y <- rowFromY(pr_rast, pr_donana$y)
  
  colnames(pr_donana)[-c(1:2)] <- paste(year(time(pr_rast)), sprintf("%02d", month(time(pr_rast))), sep = "_")
  
  
  rc_df <- as.data.frame(template, xy = T, cells = T)
  a <- rowColFromCell(template, rc_df$cell)
  rc_df$col <- a[,1]
  rc_df$row <- a[,2]
  rc_df$lyr.1 <- NULL
  rc_df$cell <- NULL
  
  # ERA5-Land 
  calculate_BM_cDM(tas_df = temp_donana,
                   pr_df = pr_donana,
                   coor_df = rc_df,
                   year_min = year_min,
                   year_max = year_max,
                   result_dir = Dir.Out,
                   temperature_type = "ERA5",
                   useRC = T)  # or FALSE for coordinates
  
}





#------------------------------------------------------------------------------------------------------------
# calculation of breeding months and consecutive dry months
#------------------------------------------------------------------------------------------------------------

# Function copied from the geosphere package - returns daylength in hours
daylength <- function(lat, doy) {
  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * 
                                                 tan(0.0086 * (doy - 186)))))
  a <- (sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(P))/(cos(lat * 
                                                                  pi/180) * cos(P))
  a <- pmin(pmax(a, -1), 1)
  DL <- 24 - (24/pi) * acos(a)
  return(DL)
}

# calculate_BM_cDM <- function(tas_df,
#                              pr_df,
#                              coor_df, 
#                              year_min = 2006, 
#                              year_max = 2025,
#                              result_dir,
#                              temperature_type = c("chelsafut", "ERA5"),
#                              useRC = T) {
#   cat('reminder: Unit of tas files should be Celcius', fill = T)
#   cat('reminder: Unit of pr files should be mm')
#   
#   cat('both file types are expected to have YYYY_MM in the file names so dates can be matched')
#   
#   if (!dir.exists(result_dir)) {dir.create(result_dir)}
#   
#   tas_time <- colnames(tas_df)[-c(1:2)]
#   pr_time <- colnames(pr_df)[-c(1:2)]
#   
#   # Check for mismatches
#   tas_time[which(!(tas_time %in% pr_time))]
#   pr_time[which(!(pr_time %in% tas_time))]
#   
#   # Define breeding probability function
#   # There are some typo's in Tablado 2012 - make sure to look at Table 2 in Tablado 2009 when looking for details
#   P_B <- function(Temp, D, delta, W) {
#     1 / (1 + exp(-(-4.542 + (0.605 * Temp) + (-0.029 * Temp^2) + (0.006 * D) + (0.017 * delta) + W)))
#   }
#   
#   # Daylength calculations
#   coor <- cbind(median(coor_df$x), median(coor_df$y))
#   pts <- vect(coor, crs = "EPSG:3035")
#   pts_wgs84 <- project(pts, "EPSG:4326")
#   coor <- as.data.frame(geom(pts_wgs84))
#   
#   dl <- daylength(coor$y, 1:365) * 60 # function returns in hours, but P_B requires daylenght in minutes
#   dL <- tapply(dl, rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31)), mean)
#   
#   # Calculate delta (day length differences)
#   ddL <- data.frame(diff_1 = dL[1] - dL[12],
#                     diff_2 = dL[2] - dL[1],
#                     diff_3 = dL[3] - dL[2],
#                     diff_4 = dL[4] - dL[3],
#                     diff_5 = dL[5] - dL[4],
#                     diff_6 = dL[6] - dL[5],
#                     diff_7 = dL[7] - dL[6],
#                     diff_8 = dL[8] - dL[7],
#                     diff_9 = dL[9] - dL[8],
#                     diff_10 = dL[10] - dL[9],
#                     diff_11 = dL[11] - dL[10],
#                     diff_12 = dL[12] - dL[11])
#   
#   # Load initial data for wet/dry calculations (November and December of previous year)
#   tas_11 <- tas_df[, c("x", "y", grep(paste0((year_min - 1), "_11"), colnames(tas_df), value = T))]
#   pr_11 <- pr_df[, c("x", "y", grep(paste0((year_min - 1), "_11"), colnames(pr_df), value = T))]
#   colnames(tas_11) <- c("x", "y", "tas")
#   colnames(pr_11) <- c("x", "y", "pr")
#   
#   suppressMessages({
#     df_11 <- left_join(
#       coor_df %>% mutate(across(c(x, y), \(v) round(v, 1))),
#       tas_11  %>% mutate(across(c(x, y), \(v) round(v, 1))),
#       by = c("x", "y")
#     ) %>% left_join(.,
#                     pr_11  %>% mutate(across(c(x, y), \(v) round(v, 1))),
#                     by = c("x", "y")              
#     )
#   })
#   
#   tas_12 <- tas_df[, c("x", "y", grep(paste0((year_min - 1), "_12"), colnames(tas_df), value = T))]
#   pr_12 <- pr_df[, c("x", "y", grep(paste0((year_min - 1), "_12"), colnames(pr_df), value = T))]
#   colnames(tas_12) <- c("x", "y", "tas")
#   colnames(pr_12) <- c("x", "y", "pr")
#   
#   suppressMessages({  
#     df_12 <- left_join(
#       coor_df %>% mutate(across(c(x, y), \(v) round(v, 1))),
#       tas_12  %>% mutate(across(c(x, y), \(v) round(v, 1))),
#       by = c("x", "y")
#     ) %>% left_join(.,
#                     pr_12  %>% mutate(across(c(x, y), \(v) round(v, 1))),
#                     by = c("x", "y")              
#     )
#   })
#   
#   # Calculate initial wet/dry conditions
#   dry_2 <- ifelse(df_11$pr < (2 * df_11$tas), TRUE, FALSE)
#   dry_1 <- ifelse(df_12$pr < (2 * df_12$tas), TRUE, FALSE)
#   
#   # Initialize result data frames
#   breed_df <- coor_df
#   dry_df <- coor_df
#   
#   b <- rep(0, nrow(coor_df))
#   
#   # Clean up temporary variables
#   rm(tas_11, pr_11, tas_12, pr_12, df_11, df_12)
#   gc()
#   
#   # Main calculation loop
#   for (yr in year_min:year_max) {
#     for (m in 1:12) {
#       
#       tas <- tas_df[, c("x", "y", grep(paste((yr), sprintf("%02d", m), sep = "_"), colnames(tas_df), value = T))]
#       pr <- pr_df[, c("x", "y", grep(paste((yr), sprintf("%02d", m), sep = "_"), colnames(pr_df), value = T))]
#       
#       colnames(tas) <- c("x", "y", "tas")
#       colnames(pr) <- c("x", "y", "pr")
#       
#       suppressMessages({
#         df <- left_join(
#           coor_df %>% mutate(across(c(x, y), \(v) round(v, 1))),
#           tas  %>% mutate(across(c(x, y), \(v) round(v, 1)))
#         ) %>% left_join(.,
#                         pr  %>% mutate(across(c(x, y), \(v) round(v, 1)))
#         )
#       })
#       
#       
#       # Get day length and delta for current month
#       D <- dL[m]
#       delta <- ddL[, m]
#       
#       # Calculate W parameter 
#       #    (W = 0 if previous two months were dry and W = -1.592 if precipitation was higher than twice the temperature in at least one of the previous two months))
#       # There are some typo's in Tablado 2012 - make sure to look at Table 2 in Tablado 2009 when looking for details
#       W <- ifelse(dry_1 & dry_2, 0, -1.592)
#       
#       # Calculate breeding probability
#       a <- ifelse(P_B(Temp = df$tas, D = D, delta = delta, W = W) >= 0.5, 1, 0)
#       a[which(is.na(a))] <- -1
#       
#       
#       # Update dry conditions for next iteration
#       dry_2 <- dry_1
#       dry_1 <- ifelse(df$pr < (2 * df$tas), TRUE, FALSE)
#       
#       # Update consecutive dry months counter
#       b <- b + 1
#       b[which(dry_1 == FALSE)] <- 0
#       b[which(is.na(W))] <- -1
#       
#       # Add results to data frames
#       breed_df <- cbind(breed_df, a)
#       dry_df <- cbind(dry_df, b)
#       
#       
#       gc()
#     }
#     
#     colnames(breed_df)[c(5:16)] <- c(1:12)
#     colnames(dry_df)[c(5:16)] <- c(1:12)
#     
#     if(useRC) {
#       breed_df <- breed_df %>% select(-c("x", "y"))
#       dry_df <- dry_df %>% select(-c("x", "y"))
#     } else {
#       breed_df <- breed_df %>% select(-c("col", "row"))
#       dry_df <- dry_df %>% select(-c("col", "row"))
#     }
#     
#     breeding_file <- file.path(result_dir, 
#                                paste0("breeding_months_", yr, '.txt'))
#     dry_months_file <- file.path(result_dir, 
#                                  paste0("consecutive_dry_months_", yr, '.txt'))
#     
#     # Write results
#     write.table(breed_df, breeding_file, quote = FALSE, row.names = FALSE)
#     write.table(dry_df, dry_months_file, quote = FALSE, row.names = FALSE)
#     
#     breed_df <- coor_df
#     dry_df <- coor_df
#     
#     gc()
#     
#   }
# }
