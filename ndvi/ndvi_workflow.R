library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(reticulate)
library(data.table)  
library(lubridate)
library(xgboost)     
library(zoo)         

env_dir = file.path("environmental_data/data")
output_dir = file.path("ndvi", "results", "model_info")
ndvi_path = file.path(env_dir, "ndvi.tif")
template_path = file.path(env_dir, "template_raster_500.tif")
krig_path = file.path(env_dir, "CDS")
precip_path = file.path(env_dir, "CDS", "precipitation.grib")

source(file.path("ndvi", "R", "quick_model_eval.R"))
source(file.path("ndvi", "R", "create_model_dataframe.R"))
source(file.path("ndvi", "R", "forecasting_function.R"))

if(!dir.exists(output_dir)) { dir.create(output_dir) }


#-------------------------------------------------------------------------------
# Load data  
#-------------------------------------------------------------------------------
template_rast <- rast(template_path)

rast_list <- create_model_rasters(template = template_rast, 
                                  ndvi_file = ndvi_path, 
                                  krig_dir = krig_path, 
                                  precip_file = precip_path)

model_data <- create_model_dataframe(ndvi_stack = rast_list$ndvi_stack, 
                                     ndvi_auto_max3 = rast_list$ndvi_neighbour, 
                                     temp_stack = rast_list$temp_stack, 
                                     precip_stack = rast_list$precip_stack)


#-------------------------------------------------------------------------------
# TRAIN MODEL - gradient boosted decision tree model
#-------------------------------------------------------------------------------

features <- c(
  "ndvi_lag1", "ndvi_lag3", "ndvi_lag12",
  "month_sin", "month_cos",
  "ndvi_max3_lag1",
  "temp", "precip"
)

train <- model_data[year < 2018]
test  <- model_data[year >= 2018]

dtrain <- xgb.DMatrix(data = as.matrix(train[, ..features]), 
                      label = train$max_ndvi)
dtest  <- xgb.DMatrix(data = as.matrix(test[, ..features]), 
                      label = test$max_ndvi)

model <- xgb.train(
  data = dtrain,
  nrounds = 100,
  objective = "reg:squarederror",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

#-------------------------------------------------------------------------------
# Quick evaluation of model metrics
#-------------------------------------------------------------------------------

# model_plots <- basic_eval()

#-------------------------------------------------------------------------------
# save model for easier forecasting later on
#-------------------------------------------------------------------------------

# model
xgb.save(model, file.path(output_dir, "ndvi_xgb_model.json"))

# residuals
preds <- predict(model, as.matrix(train[, ..features]))
residuals <- train$max_ndvi - preds

saveRDS(residuals, file.path(output_dir, "residuals.rds"))

# input data
writeRaster(rast_list$ndvi_stack, file.path(output_dir, "ndvi_stack.tif"), overwrite = TRUE)
writeRaster(rast_list$ndvi_neighbour, file.path(output_dir, "neighbour_stack.tif"), overwrite = TRUE)
writeRaster(rast_list$temp_stack, file.path(output_dir, "temp_stack.tif"), overwrite = TRUE)
writeRaster(rast_list$precip_stack, file.path(output_dir, "precip_stack.tif"),   overwrite = TRUE)

saveRDS(time(rast_list$ndvi_stack), file.path(output_dir, "ndvi_time.rds"))
saveRDS(time(rast_list$ndvi_neighbour), file.path(output_dir, "neighbour_time.rds"))
saveRDS(time(rast_list$temp_stack), file.path(output_dir, "temp_time.rds"))
saveRDS(time(rast_list$precip_stack), file.path(output_dir, "precip_time.rds"))


#-------------------------------------------------------------------------------
# FORECAST
#-------------------------------------------------------------------------------


ndvi_forecast <- forecast_ndvi_raster_multivar(
  model_dir = output_dir, 
  climate_type = "mean",
  stochastic = FALSE
)

ndvi_mean <- global(ndvi_forecast, fun = "mean", na.rm = TRUE)
ndvi_mean$time <- time(ndvi_forecast)



pred_df <- lapply(as.list(1:10), function(x) {
  ndvi_forecast <- forecast_ndvi_raster_multivar(
    model_dir = output_dir, 
    climate_years = c(2015, 2016),
    stochastic = TRUE
  )
  ndvi_df <- global(ndvi_forecast, fun = "mean", na.rm = TRUE)
  ndvi_df$time <- time(ndvi_forecast)
  return(ndvi_df)
}) %>% bind_rows(.id = "rep")


ggplot() +
  geom_line(data = pred_df, aes(x = time, y = mean, group = rep), colour = "green2") +
  geom_line(data = ndvi_mean,aes(x = time, y = mean), colour = "black" ) +
  theme_minimal()








