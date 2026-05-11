library(dplyr)
library(terra)

# Set paths---------------------------------------------------------------------
data_dir = file.path("environmental_data", "data")
border_path = file.path(data_dir, "Donana_DT_border", "Limite_Don╠âana.shp")
template_path = file.path(data_dir, "template_raster_500.tif")
dem_path = file.path(data_dir, "Copernicus_GLO90_Europe_250m.tif")

# load functions----------------------------------------------------------------
# this file just assigns two objects: API_User and API_Key but stored in data to keep it private
source(file.path(data_dir, "CDS_API_info.R"))
source(file.path("environmental_data", "functions", "Download_and_Krig_temperature.R"))


#-------------------------------------------------------------------------------
# Create the raster template
#-------------------------------------------------------------------------------
border_vec <- vect(border_path)
border_vec <- buffer(border_vec, 1000)

# Create template rasters with 500m resolution
DT_template <- rast(
  xmin = xmin(border_vec), 
  xmax = xmax(border_vec),
  ymin = ymin(border_vec), 
  ymax = ymax(border_vec),
  resolution = c(500, 500),
  crs = crs(border_vec)
)
values(DT_template) <- 0

writeRaster(DT_template, template_path, 
            overwrite = TRUE, NAflag = -9999)

#-------------------------------------------------------------------------------
# Download and Krig temperature data
#-------------------------------------------------------------------------------

download_and_krig_temp(template_path = template_path, dem_path = dem_path, 
                       output_dir = file.path(data_dir, "CDS"))


#-------------------------------------------------------------------------------
# Download precipitation data
#-------------------------------------------------------------------------------

cat('still need to figure out how to do this only for the months not yet downloaded')

if(!file.exists(file.path(data_dir, "CDS", "precipitation.grib"))){
  cat('downloading precipitation data - this might take a while')
  # use_python("C:/Users/z1512834z/AppData/Local/Programs/Python/Python314/python.exe")
  source_python('environmental_data/functions/download_precipitation_CDS.py')
}

#-------------------------------------------------------------------------------
# Download and clean NDVI
#-------------------------------------------------------------------------------

# More info in ndvi/python/NDVI-Download-Pipeline/README.md

# I had problems with a conflicting proj db. Creating a clean pyton environment first to avoid conflict
env_dir <- "python_env"
python_exe <- file.path(env_dir, "Scripts", "python.exe")

if (!file.exists(python_exe)) {
  system(sprintf('python -m venv "%s"', env_dir))

  system(sprintf('"%s" -m pip install --upgrade pip', python_exe))
  # system(sprintf('"%s" -m pip install rasterio pyproj', python_exe))

  # Installing required modules
  system2(python_exe, args = c("-m", "pip", "install", "-r", "environmental_data/functions/NDVIDownloadPipeline/requirements.txt"))

}
# 
# start <- ymd("2005-01-01")
# end   <- floor_date(Sys.Date(), "month") - days(1)
# seq(start, end, by = "month")
# 

system2(
  python_exe,
  args = c("environmental_data/functions/NDVIDownloadPipeline/pipeline.py",
           "--output-dir", "environmental_data/data/NDVI/",
           "--shape-path", "environmental_data/data/Donana_DT_border/Limite_Don╠âana.shp",
           "--date", "2004-01-01/2026-04-30")
)


n <- rast("environmental_data/functions/NDVIDownloadPipeline/output/ndvi.nc")



#-------------------------------------------------------------------------------
# Flooding data
#-------------------------------------------------------------------------------

# Download the flooding maps on the WCS ----------------------------------------
system("environmental_data/functions/download_flooding_wcs.py")

# The most recent files are not on there yet, so download ----------------------

fld_files <- list.files(file.path(data_dir, "LAST_Inundacion"), full.names = T)

dates <- stringr::str_extract(fld_files, "\\d{8}")
## Below round the date of taking floods that happen throughout the month, the flood is taken into
## account at the load of the closest "new month"
dates <- as.Date(dates, format = "%Y%m%d")
yrs <- unique(year(dates))
yrs <- yrs[which(yrs > 1999)]
mnts <- c(1:12)

for(y in yrs){
  for(m in mnts){
    try({
      
      r_files <- fld_files[which(year(dates) == y & month(dates) == m)]
      
      if (length(r_files) == 0) next
      
      r_list <- lapply(r_files, function(f) {
        r <- rast(f)
        r <- project(r, crs(DT_template))
        r <- resample(r, DT_template, method = "max")
        return(r)
      })
      
      r_max <- round(max(rast(r_list), na.rm = TRUE), digits = 0) 
      
      writeRaster(r_max, 
                  file.path(data_dir, "LAST_monthly_floodmaps", paste0(y, "_", sprintf("%02d", m), ".asc")),
                  overwrite = T, NAflag = 0)
    })
  }
}



#-------------------------------------------------------------------------------
# Visualise some data
#-------------------------------------------------------------------------------

source(file.path("environmental_data", "functions", "dry_wet_years.R"))

temp_files = list.files(file.path(data_dir, "CDS"), 
                        pattern = "Kriged.nc$", full.names = T)
temp_stack = lapply(temp_files, rast) %>% rast(.)
temp_stack = project(temp_stack, crs(DT_template))
temp_stack = resample(temp_stack, DT_template, method = "bilinear")
temp_stack = temp_stack - 273.15

# Precipitation -------------------------------------------------------------------

precip_stack = rast(file.path(data_dir, "CDS", "precipitation.grib"))
precip_stack = project(precip_stack, crs(DT_template))
precip_stack = resample(precip_stack, DT_template, method = "bilinear")

# transform precipitation from mean m/day to total mm/month
dates <- time(precip_stack)
scaling <- 1000 * days_in_month(dates)
precip_scaled <- precip_stack

for (i in 1:nlyr(precip_stack)) {
  precip_scaled[[i]] <- precip_stack[[i]] * scaling[i]
}



years <- dry_wet_years(r_temp = temp_stack,
                       r_precip = precip_scaled)

year_df <- years %>% group_by(year) %>% 
  summarise(z_temp = mean(z_temp, na.rm = T),
            z_precip = mean(z_precip, na.rm = T)) %>%
  filter(year != 2004)



year_plot <- ggplot(year_df, aes(x = z_temp, y = z_precip, label = year)) +
  geom_text(size = 3) +
  annotate("segment", 
           x = -1, xend = 1,
           y = 0, yend = 0,
           arrow = arrow(length = unit(0.25, "cm"),
                         ends = "both", type = "closed"),
           color = "black") +
  annotate("text", x = 1, y = -0.2,
           label = "Hot", vjust = -0.5,
           color = "red") +
  annotate("text", x = -1, y = -0.2,
           label = "Cold", vjust = -0.5,
           color = "orange") +
  annotate("segment",
           x = 0, xend = 0,
           y = -0.7, yend = 1,
           arrow = arrow(length = unit(0.25, "cm"),
                         ends = "both", type = "closed"),
           color = "black")  +
  annotate("text", x = 0.2, y = 1,
           label = "Wet", vjust = -0.5,
           color = "blue") +
  annotate("text", x = -0.2, y = -0.7,
           label = "Dry", vjust = -0.5,
           color = "orange") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) +
  coord_cartesian()

ggsave(year_plot, filename = file.path(data_dir, "wet_dry_year_distribution.png"),
       width = 5, height = 5)




