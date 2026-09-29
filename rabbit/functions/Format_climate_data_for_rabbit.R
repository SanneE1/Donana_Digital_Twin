
library(dplyr)

##------------------------------------------------------------------------------------------------------------
## Use downloaded data to calculate breeding months and consecutive dry months
##------------------------------------------------------------------------------------------------------------

# Function copied from the geosphere package
daylength <- function(lat, doy) {
  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * 
                                                 tan(0.0086 * (doy - 186)))))
  a <- (sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(P))/(cos(lat * 
                                                                  pi/180) * cos(P))
  a <- pmin(pmax(a, -1), 1)
  DL <- 24 - (24/pi) * acos(a)
  return(DL)
}

calculate_BM_cDM <- function(temp,
                             precip,
                             K = T,
                             coord_rast, 
                             year_min = 2002, 
                             year_max = 2025,
                             result_dir) {
  
  if (!dir.exists(result_dir)) {dir.create(result_dir)}
  
  if(K){temp <- temp - 273.15}
  
  tas_time <- floor_date(time(temp), "month")
  pr_time  <- floor_date(time(precip), "month")
  
  # Check for mismatches
  tas_time[which(!(tas_time %in% pr_time))]
  pr_time[which(!(pr_time %in% tas_time))]
  
  # Define breeding probability function
  P_B <- function(Temp, D, delta, W) {
    1 / (1 + exp(-(-4.542 + (0.605 * Temp) + (-0.029 * Temp^2) + (0.006 * D) + (0.017 * delta) + W)))
  }
  
  # Daylength calculations
  mL <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mE <- cumsum(mL)
  
  dL <- lapply(as.list(c(1:12)), function(s) {
    crds(coord_rast, df=T) %>%
      rowwise() %>%
      mutate(dayL = mean(daylength(y, c((mE[s] - mL[s]):mE[s])) * 60),
             .keep = "none")
  }) %>% bind_cols()
  
  # Calculate delta (day length differences)
  ddL <- data.frame(diff_1 = dL[,1] - dL[,12],
                    diff_2 = dL[,2] - dL[,1],
                    diff_3 = dL[,3] - dL[,2],
                    diff_4 = dL[,4] - dL[,3],
                    diff_5 = dL[,5] - dL[,4],
                    diff_6 = dL[,6] - dL[,5],
                    diff_7 = dL[,7] - dL[,6],
                    diff_8 = dL[,8] - dL[,7],
                    diff_9 = dL[,9] - dL[,8],
                    diff_10 = dL[,10] - dL[,9],
                    diff_11 = dL[,11] - dL[,10],
                    diff_12 = dL[,12] - dL[,11])
  
  # Load initial data for wet/dry calculations (November and December of previous year)
  tas_11 <- subset(temp, tas_time == as.Date(paste(year_min-1, "11", "01", sep = "-"))) 
  tas_11 <- project(tas_11, crs(coord_rast), method = "bilinear")
  tas_11 <- terra::extract(tas_11, crds(coord_rast))
  
  pr_11 <- subset(precip, pr_time == as.Date(paste(year_min-1, "11", "01", sep = "-"))) 
  pr_11 <- project(pr_11, crs(coord_rast), method = "bilinear")
  pr_11 <- terra::extract(pr_11, crds(coord_rast)) * 30 * 1000
  
  tas_12 <- subset(temp, tas_time == as.Date(paste(year_min-1, "12", "01", sep = "-"))) 
  tas_12 <- project(tas_12, crs(coord_rast), method = "bilinear")
  tas_12 <- terra::extract(tas_12, crds(coord_rast))
  
  pr_12 <- subset(precip, pr_time == as.Date(paste(year_min-1, "12", "01", sep = "-"))) 
  pr_12 <- project(pr_12, crs(coord_rast), method = "bilinear")
  pr_12 <- terra::extract(pr_12, crds(coord_rast)) * 31 * 1000
  
  # Calculate initial wet/dry conditions
  dry_2 <- ifelse(pr_11 < (2 * tas_11), TRUE, FALSE)
  dry_1 <- ifelse(pr_12 < (2 * tas_12), TRUE, FALSE)
  
  # Initialize result data frames
  coords_df <- crds(coord_rast, df=T)
  colRow_df <- data.frame(col = colFromCell(coord_rast, cells(coord_rast)),
                          row = rowFromCell(coord_rast, cells(coord_rast)))
  breed_df_empty <- colRow_df
  dry_df_empty <- colRow_df
  
  breed_df <- breed_df_empty
  dry_df <- dry_df_empty
  
  b <- rep(0, length(dry_df[, 1]))
  
  # Clean up temporary variables
  rm(tas_11, pr_11, tas_12, pr_12)
  gc()
  
  # Main calculation loop
  for (yr in year_min:year_max) {
    for (m in 1:12) {
      
      tas <- subset(temp, tas_time == as.Date(paste(yr, sprintf("%02d", m), "01", sep = "-"))) 
      tas <- project(tas, crs(coord_rast), method = "bilinear")
      tas <- terra::extract(tas, crds(coord_rast))
      
      pr <- subset(precip, pr_time == as.Date(paste(yr, sprintf("%02d", m), "01", sep = "-"))) 
      pr <- project(pr, crs(coord_rast), method = "bilinear")
      pr <- terra::extract(pr, crds(coord_rast)) * days_in_month(m) * 1000
      
      
      # Get day length and delta for current month
      D <- dL[, m]
      delta <- ddL[, m]
      
      # Calculate W parameter
      W <- ifelse(dry_1 & dry_2, 0, -1.592)
      
      # Calculate breeding probability
      a <- ifelse(P_B(Temp = tas, D = D, delta = delta, W = W) >= 0.5, 1, 0)
      a[which(is.na(a))] <- -1
      
      # Update dry conditions for next iteration
      dry_2 <- dry_1
      dry_1 <- ifelse(pr < (2 * tas), TRUE, FALSE)
      
      # Update consecutive dry months counter
      b <- b + 1
      b[which(dry_1 == FALSE)] <- 0
      b[which(is.na(W))] <- b-1
      
      # Add results to data frames
      breed_df <- cbind(breed_df, a)
      dry_df <- cbind(dry_df, b)
      
      
      gc()
    }
    
    colnames(breed_df)[c(3:14)] <- c(1:12)
    colnames(dry_df)[c(3:14)] <- c(1:12)
    
    which(is.na(breed_df))
    
    breeding_file <- file.path(result_dir, 
                               paste0("breeding_months_", yr, '.txt'))
    dry_months_file <- file.path(result_dir, 
                                 paste0("consecutive_dry_months_", yr, '.txt'))
    
    # Write results
    write.table(breed_df, breeding_file, quote = FALSE, row.names = FALSE)
    write.table(dry_df, dry_months_file, quote = FALSE, row.names = FALSE)
    
    breed_df <- breed_df_empty
    dry_df <- dry_df_empty
    
    gc()
    
  }
  
  
  
}
