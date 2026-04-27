
library(terra)
library(tidyverse)
library(tidyterra)
library(gganimate)

source("R/Rasterize_output_maps.R")

DT_border <- vect("data/original_data/Donana_DT_border/Limite_Don╠âana.shp")
hab_map <- rast("data/GIS_maps/Rabbit_HabitatMap_500_Donana_Fordham_2013.asc")
DT_border <- project(DT_border, crs(hab_map))

if(exists("run_files")) {
  map_files <- list.files(run_files, full.names = T)
} else {map_files <- list.files("results/forecast_summary", full.names = T)}

maps <- lapply(map_files, function(x) csvToRaster(x, hab_map))
tif_names <- file.path("results", "GIS_summary_maps", 
                       str_replace_all(basename(map_files), pattern = ".csv", replacement = ".tif"))

for(i in c(1:length(tif_names))){
  writeRaster(maps[[i]], filename = tif_names[[i]], overwrite = TRUE)
}

map_list <- rast(maps)
names(map_list) <- sapply(map_files, function(x) {
  
  a <- str_split(x, "[[:punct:]]")[[1]][c(7,8,9)]
  a[2] <- sprintf("%02d", as.integer(a[2])) 
  paste(a, collapse = "-")
})

mean_maps <- map_list[[sort(names(map_list)[grepl('mean', names(map_list))])]]

plots <- lapply(as.list(names(mean_maps)), function(i) {
  ggplot() +
    geom_spatraster(data = mean_maps[[i]]) +
    geom_spatvector(data = DT_border, aes(fill = NA), colour = "black", size = 1) + 
    scale_fill_viridis_c(
      trans = "pseudo_log",
      limits = c(0,ceiling(max(minmax(mean_maps[[names(mean_maps)[grepl('mean', names(mean_maps))]]])))), na.value = NA) +
    labs(title = paste("file:", i )) + theme_classic()
})

names(plots) <- names(mean_maps)


## Create a gif of mean plots
library(magick)
dir.create("results/frames", showWarnings = FALSE)

for (i in seq_along(plots)) {
  ggsave(
    filename = sprintf("results/frames/frame_%02d.png", i),
    plot = plots[[i]],
    width = 6,
    height = 4,
    dpi = 150
  )
}

gif <- image_read(list.files("results/frames", full.names = TRUE)) |>
  image_animate(fps = 1)  # 1 second per plot

image_write(gif, "plots.gif")


# Calibration comparison
transects <- vect("data/original_data/Rabbit_donana_KAI_PacoCarro/Transect_oryctolagus.kml")
vect_length <- terra::perim(transects)/1000  # Perim gives m, but need km for KAI - not that it matters as we're working with relative index

area_trans <- data.frame(ID = c(1:37),
                         transect = c("Abalario", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Hinojos",
                                      "Hinojos-Guadiamar", "Marismillas", "Muro", "Matochal", "Puntal",
                                      "RBD-este", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea",  
                                      "Sabinar-Mogea", "RBD-este", "RBD-este", "Algaida-Sotos", "Algaida-Sotos", 
                                      "Algaida-Sotos", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Coto del Rey", 
                                      "Coto del Rey", "Coto del Rey", "Coto del Rey", "Coto del Rey", "Hinojos", 
                                      "Hinojos", "Marismillas", "Marismillas", "Marismillas", "Puntal",
                                      "Puntal", "Puntal"
                         ))

length_trans <- cbind(area_trans, vect_length) %>% 
  group_by(transect) %>%
  summarise(trans_length = sum(vect_length))

transects <- project(transects, crs(hab_map))

obs_data <- read.csv("data/original_data/Rabbit_donana_KAI_PacoCarro/KAI_Rabbit_Night_2024_v1.csv") %>% 
  mutate(date = as.Date(Fecha, "%d/%m/%Y")) %>%
  select(-c(Fecha, Fecha.1, Media.END)) %>%
  pivot_longer(cols = -date,
               names_to = "transect", 
               values_to = "obs_KAI") %>%
  mutate(transect = gsub("\\.", " ", transect))


sim_transect <- terra::extract(map_list, transects) %>%
  left_join(., area_trans) %>%
  select("transect", contains("mean"), contains("lower"), contains("upper")) %>%
  group_by(transect) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = T))) %>%
  pivot_longer(
    cols = -transect,
    names_to = c("date", ".value"),
    names_sep = "-(?=mean|lower|upper)"
  ) %>% 
  mutate(date = lubridate::ym(date),
         transect = gsub("\\-", " ", transect)) %>%
  arrange(date)


ggplot() + 
  geom_ribbon(data = sim_transect, aes(x = date, ymin = lower, ymax = upper), colour = "green") + 
  geom_line(data = sim_transect, aes(x = date, y = mean), colour = "darkgreen") +
  geom_line(data = obs_data, aes(x = date, y = obs_KAI*1000), colour = "red") +
  facet_wrap(vars(transect), scales = "free") + 
  scale_y_continuous(name = "Simulated abundance", 
                     sec.axis = sec_axis(~./1000, name = "Observed KAI")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.left = element_text(color = "darkgreen"),
        axis.text.y.left = element_text(color = "darkgreen"))


