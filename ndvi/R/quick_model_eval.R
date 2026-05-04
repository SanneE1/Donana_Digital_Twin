library(pdp)


basic_eval <- function() {
  test[, pred := predict(model, as.matrix(test[, ..features]))]
  test[, time := as.Date(paste(test[,year], test[,month], "15", sep = "/"))]
  test[, residual := max_ndvi - pred]
  
  rmse <- sqrt(mean((test$max_ndvi - test$pred)^2))
  mae  <- mean(abs(test$max_ndvi - test$pred))
  r2   <- cor(test$max_ndvi, test$pred)^2
  
  cat("RMSE:", rmse, "\n")
  cat("MAE :", mae, "\n")
  cat("R²  :", r2, "\n")
  
  
  pred_rast <- rast(test[, .(x, y, time, pred)], type = "xylz")
  
  eval1 <- ggplot(test, aes(x = residual)) +
    geom_histogram(bins = 50) +
    theme_minimal() + ggtitle("Residual distribution")
  
  eval2 <- ggplot(test, aes(x = pred, y = residual)) +
    geom_point(alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red") +
    theme_minimal() + ggtitle("Predicted vs Residual")
  
  
  test_monthly <- test[, .(
    ndvi = mean(max_ndvi),
    pred = mean(pred)
  ), by = time]
  
  eval3 <- ggplot(test_monthly, aes(x = time)) +
    geom_line(aes(y = ndvi, color = "observed")) +
    geom_line(aes(y = pred, color = "predicted")) +
    theme_minimal() +
    scale_color_manual(name = "data",
                       values = c("red", "blue")) +
    labs(title = "Temporal dynamics: observed vs predicted")
  
  spatial_error <- test[, .(
    rmse = sqrt(mean((max_ndvi - pred)^2))
  ), by = .(x, y)]
  
  r_err <- rast(spatial_error[, .(x, y, rmse)], type = "xyz")
  
  eval4 <- ggplot() +
    geom_spatraster(data = r_err) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = "Spatial distribution of rmse")
  
  plot1 <- eval1 + eval2 + eval3 + eval4 + plot_layout(ncol = 2)
  
  
  importance <- xgb.importance(model = model)
  print(importance)
  
  xgb.plot.importance(importance)
  
  
  coeff_list <- list(
    "ndvi_lag1", "ndvi_lag3", "ndvi_lag12",
    "ndvi_max3_lag1", "temp", "precip"
  )
  
  parts_list<- lapply(coeff_list, function(x) {
    p_part <- partial(
      model,
      pred.var = x,
      train = as.data.frame(train[, ..features])
    )
  })
  
  
  plot2 <- (ggplot(data = parts_list[[1]]) +
              geom_line(aes(x = ndvi_lag1, y = yhat), color = "darkred", linewidth = 1) +
              ylab("Predicted NDVI") + 
              theme_minimal()) +
    (ggplot(data = parts_list[[2]]) +
       geom_line(aes(x = ndvi_lag3, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal()) +
    (ggplot(data = parts_list[[3]]) +
       geom_line(aes(x = ndvi_lag12, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal()) +
    (ggplot(data = parts_list[[4]]) +
       geom_line(aes(x = ndvi_max3, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal()) +
    (ggplot(data = parts_list[[5]]) +
       geom_line(aes(x = ndvi_max3_lag1, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal()) +
    (ggplot(data = parts_list[[6]]) +
       geom_line(aes(x = temp, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal()) +
    (ggplot(data = parts_list[[7]]) +
       geom_line(aes(x = precip, y = yhat), color = "darkred", linewidth = 1) +
       ylab("Predicted NDVI") + 
       theme_minimal())  +
    plot_layout(ncol = 2)
  
  return(list(plot_residuals = plot1,
              plot_partial = plot2))
  
}