#
# Lynx observations and model input files for the Doñana lynx components, from the PVA project
# (C:/Users/z1512834z/Documents/01_Natura_Connect/04_Lynx_PVA, not in this repository)
#
#   data/lynx/census/donana_breeding_females.csv
#       year, breeding_females, source. 2002-2020: IUCN Doñana sizes (mature individuals) / 2,
#       which matches LynxConnect where both exist (2021: IUCN 46, LynxConnect 23 breeding
#       females); 2021-2024: LynxConnect census, Doñana-Aljarafe, "Hembras reproductoras".
#   data/lynx/presence/presence_<year>.gpkg
#       annual lynx presence (PVA presence_vectors, 2002-2018 dissolved 1 km areas, 2020-2022 1 km
#       cells; no 2019), cropped to the large template (+ 5 km)
#   data/lynx/Lynx_demography_values_donana.txt
#       the PVA rates (Lynx_demography_values_correctedIPM.txt: min_rep_age 2, BOM removed);
#       survival, Tsize, rep_prob and surv_disp_rho are calibrated
#
# Run from the Donana_Digital_Twin project root.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

pva <- "C:/Users/z1512834z/Documents/01_Natura_Connect/04_Lynx_PVA/data"
out <- file.path("data", "lynx")
dir.create(file.path(out, "census"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "presence"), showWarnings = FALSE)
template <- rast(file.path("data", "environmental_data", "template_raster_500_large.tif"))

iucn <- read.csv(file.path(pva, "original_data", "Population_sizes_IUCN.csv"), check.names = FALSE,
                 fileEncoding = "UTF-8-BOM")
lc <- read.csv(file.path(pva, "original_data", "2025.08.06_LynxConnectWebsiteCensusNumber.csv"),
               fileEncoding = "latin1", check.names = FALSE)
lc <- lc[grepl("Aljarafe", lc[[2]]), ]
counts <- rbind(
  data.frame(year = iucn$Year, breeding_females = iucn[["Doñana"]] / 2, source = "IUCN Doñana mature individuals / 2"),
  data.frame(year = lc$Year, breeding_females = lc[["Hembras reproductoras"]], source = "LynxConnect Doñana-Aljarafe"))
counts <- counts[!(counts$source != "LynxConnect Doñana-Aljarafe" & counts$year %in% lc$Year), ]
counts <- counts[order(counts$year), ]
write.csv(counts, file.path(out, "census", "donana_breeding_females.csv"), row.names = FALSE)
print(counts)

area <- buffer(as.polygons(ext(template), crs = crs(template)), 5000)
for (f in list.files(file.path(pva, "GIS_maps", "presence_vectors"), "\\.shp$", full.names = TRUE)) {
  v <- crop(project(vect(f), crs(template)), area)
  if (nrow(v) == 0) next
  v <- aggregate(v)            # one (multi)polygon: only the area matters
  writeVector(v, file.path(out, "presence", sprintf("presence_%s.gpkg", sub("\\.shp$", "", basename(f)))),
              overwrite = TRUE)
}
cat("presence years:", sub("presence_(\\d+).*", "\\1", list.files(file.path(out, "presence"))), "\n")

rates <- readLines(file.path(pva, "model_input", "Lynx_demography_values_correctedIPM.txt"), encoding = "UTF-8")
rates[1] <- sub("^\ufeff", "", rates[1])
# max_age 15: the PVA survival is 0 from age 16 on. write_lynx_rates() sets the calibrated adult
# survival for ages 2..max_age and 0 above, so no lynx gets older than the survival table (ages
# 0-21; age 22 crashed the IBM)
rates <- trimws(rates)
rates[startsWith(rates, "max_age ")] <- "max_age 15"
writeLines(rates, file.path(out, "Lynx_demography_values_donana.txt"))
