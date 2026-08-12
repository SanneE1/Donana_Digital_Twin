library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(reticulate)
library(data.table)  
library(lubridate)
library(xgboost)     
library(zoo)         


# Settings
test_data_prop = 0.2

env_dir = file.path("data", "environmental_data")
output_dir = file.path("ndvi", "results")

model_dir = file.path(output_dir, "model_info")
pred_dir = file.path(output_dir, "predictions")

ndvi_path = file.path(env_dir, "ndvi.tif")
template_path = file.path(env_dir, "template_raster_500.tif")
krig_path = file.path(env_dir, "CDS")
precip_path = file.path(env_dir, "CDS", "precipitation.nc")

sapply(list.files(file.path("ndvi", "functions"), full.names = T), source)

if(!dir.exists(output_dir)) { dir.create(output_dir, recursive = T) }
if(!dir.exists(pred_dir)) { dir.create(pred_dir, recursive = T) }
if(!dir.exists(model_dir)) { dir.create(model_dir, recursive = T) }

#-------------------------------------------------------------------------------
# Load data  
#-------------------------------------------------------------------------------
template_rast <- rast(template_path)
ndvi_files <- list.files(ndvi_path)

rast_list <- create_model_rasters(template = template_rast, 
                                  ndvi_file = ndvi_path, 
                                  krig_dir = krig_path, 
                                  precip_file = precip_path)

model_data <- create_model_dataframe(ndvi_stack = rast_list$ndvi_stack, 
                                     ndvi_auto_max3 = rast_list$ndvi_neighbour, 
                                     temp_stack = rast_list$temp_stack, 
                                     precip_stack = rast_list$precip_stack)

n_years_test <- round((max(model_data$year) - min(model_data$year)) * test_data_prop)
cut_off_year <- max(model_data$year) - n_years_test

cat('\nUsing the last ', n_years_test, 'years for the test data, everything before', 
    cut_off_year, 'will be used to train the model')

#-------------------------------------------------------------------------------
# TRAIN MODEL - gradient boosted decision tree model
#-------------------------------------------------------------------------------

features <- c(
  "ndvi_lag1", "ndvi_lag3", "ndvi_lag12",
  "month_sin", "month_cos",
  "ndvi_max3_lag1",
  "temp", "precip"
)

train <- model_data[year < cut_off_year]
test  <- model_data[year >= cut_off_year]

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

model_plots <- basic_eval()

#-------------------------------------------------------------------------------
# save model for forecasting later on
#-------------------------------------------------------------------------------

# model
xgb.save(model, file.path(model_dir, "ndvi_xgb_model.json"))

# residuals
preds <- predict(model, as.matrix(train[, ..features]))
residuals <- train$max_ndvi - preds

saveRDS(residuals, file.path(model_dir, "residuals.rds"))

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
# Save climate objects for forecasts
#-------------------------------------------------------------------------------

temp_summary <- list(mean = tapp(rast_list$temp_stack, "month", "mean"),
                     sd = tapp(rast_list$temp_stack, "month", "sd"))


precip_summary <- list(mean = tapp(rast_list$precip_stack, "month", "mean"),
                       sd = tapp(rast_list$precip_stack, "month", "sd"))


temp_rast <- list(warm = (2 * temp_summary$sd) + temp_summary$mean,
                  mean = temp_summary$mean,
                  cold = (-2 * temp_summary$sd) + temp_summary$mean)

precip_rast <- list(wet = (2 * precip_summary$sd) + precip_summary$mean,
                    mean = precip_summary$mean,
                    dry = (-2 * precip_summary$sd) + precip_summary$mean)

writeRaster(temp_rast$warm, file.path(output_dir, "temp_warm_values.tif"), overwrite = TRUE)
writeRaster(temp_rast$mean, file.path(output_dir, "temp_mean_values.tif"), overwrite = TRUE)
writeRaster(temp_rast$cold, file.path(output_dir, "temp_cold_values.tif"), overwrite = TRUE)
writeRaster(precip_rast$wet, file.path(output_dir, "precip_wet_values.tif"),   overwrite = TRUE)
writeRaster(precip_rast$mean, file.path(output_dir, "precip_mean_values.tif"), overwrite = TRUE)
writeRaster(precip_rast$dry, file.path(output_dir, "precip_dry_values.tif"),   overwrite = TRUE)

saveRDS(time(temp_rast$warm), file.path(output_dir, "temp_warm_time.rds"))
saveRDS(time(temp_rast$mean), file.path(output_dir, "temp_mean_time.rds"))
saveRDS(time(temp_rast$cold), file.path(output_dir, "temp_cold_time.rds"))
saveRDS(time(precip_rast$wet), file.path(output_dir, "precip_wet_time.rds"))
saveRDS(time(precip_rast$mean), file.path(output_dir, "precip_mean_time.rds"))
saveRDS(time(precip_rast$dry), file.path(output_dir, "precip_dry_time.rds"))
