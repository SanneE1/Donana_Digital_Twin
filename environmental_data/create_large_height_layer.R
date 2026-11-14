#
# Static vegetation height (m) on the large (lynx) template, for the rabbit K map outside the park
#
#   vegetation_height_500_large.tif   mean height per 500 m cell
#
# Inside the LiDAR footprint (vegetation_height/merged_veg_height_donana_NA_interpolated.tif, 2.5 m)
# the mean LiDAR height, as the rabbit calibration uses. Elsewhere the ETH Global Canopy Height 2020
# (10 m, exported at 100 m by backcast/ndvi_pipeline/height_water_gee.py height), rescaled to the
# LiDAR with a linear fit over the cells that have both. ETH is about 2.7 times higher than the
# LiDAR in the park (mean 2.1 vs 0.8 m) but ranks the cells the same way (500 m: r 0.90, R2 0.81).
#
# Run from the Donana_Digital_Twin project root after create_large_template.R.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

data_dir <- file.path("data", "environmental_data")
template <- rast(file.path(data_dir, "template_raster_500_large.tif"))[[1]]

lidar <- rast(file.path(data_dir, "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif"))
eth   <- rast(file.path(data_dir, "height_large", "height100_eth.tif"))
eth[eth == -9999] <- NA
crs(eth) <- crs(template)   # same EPSG:25829, the GEE export writes a different WKT

lidar500 <- resample(lidar, template, method = "average")
eth500   <- resample(eth, template, method = "average")
# Cells only partly covered by the LiDAR are left to ETH: share of 2.5 m LiDAR cells per 500 m cell
cover <- resample(!is.na(lidar), template, method = "average")
lidar500[is.na(cover) | cover < 0.9] <- NA

d <- na.omit(data.frame(lidar = values(lidar500, mat = FALSE), eth = values(eth500, mat = FALSE)))
fit <- lm(lidar ~ eth, d)
cat(sprintf("LiDAR = %.3f + %.3f * ETH (n %d cells, R2 %.2f)\n", coef(fit)[1], coef(fit)[2], nrow(d), summary(fit)$r.squared))

height <- cover(lidar500, coef(fit)[1] + coef(fit)[2] * eth500)
# ETH has no value over water (30 cells, all salines/water, K = 0 in rabbit_no_habitat_500_large.tif):
# bare ground (the fit at ETH = 0)
n_na <- sum(!is.na(values(template, mat = FALSE)) & is.na(values(height, mat = FALSE)))
cat(n_na, "land cells without LiDAR or ETH set to", round(coef(fit)[1], 3), "m\n")
height <- mask(cover(height, init(template, coef(fit)[1])), template)
names(height) <- "vegetation_height"
writeRaster(height, file.path(data_dir, "vegetation_height_500_large.tif"), overwrite = TRUE)
cat(sprintf("cells from LiDAR %d, from ETH %d\n", sum(!is.na(values(mask(lidar500, template), mat = FALSE))),
            sum(!is.na(values(height, mat = FALSE))) - sum(!is.na(values(mask(lidar500, template), mat = FALSE)))))
