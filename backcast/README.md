# Rabbit backcast

Reconstructs monthly rabbit densities as far back as the environmental data allow, with the
calibrated rabbit model (`rabbit/results/posterior_samples.csv`, calibrated on 2017–2024).
Builds on the forecast workflow's corrected inputs, which it only reads:
`data/environmental_data/CDS/kriged_temperature_v2` and `data/environmental_data/monthly_ndvi_v2`
(run `environmental_data/rebuild_monthly_ndvi.R` first). Everything it writes goes to
`data/backcast` and `rabbit/results/backcast_<start year>`.

NB: the posterior was calibrated on the earlier inputs (temperature labelled one month early, NDVI
with gaps). Until the model is recalibrated on the corrected inputs, backcast results are
provisional.

| Phase | Simulated | NDVI | Run |
|---|---|---|---|
| 1 | 2015 → last full year (2015 = spin-up) | Sentinel-2 (from July 2015) | `Workflow_backcast_data.R` with `phase <- 1`, then `rabbit_backcast(start_year = 2015)` |
| 2 | 1985 → last full year (1985 = spin-up) | Landsat 5/7/8 adjusted to Sentinel-2 before July 2015 | `Workflow_backcast_data.R` with `phase <- 2`, `harmonise_landsat_to_s2.R`, then `rabbit_backcast(start_year = 1985, ndvi_mean_file = "data/backcast/monthly_ndvi_harmonised/mean.tif", ndvi_sd_file = "data/backcast/monthly_ndvi_harmonised/sd.tif")` |

Run everything from the project root; `rabbit_backcast()` is defined in `Workflow_backcast_rabbit.R`.

## Backcast from 1985 with the rabbit components (large template)

`Workflow_backcast_components.R` prepares the inputs and then runs the two LifeWatch components in
`../Components`: `donana-rabbit-calibration` on the park (cut out of the large inputs), simulated
from 1985 so that the RHD arrival dates (1990, 2012) take effect and fitted to the census
2005-2024, and `donana-rabbit-forecasting` with the observed NDVI and climate
(`useNdviForecast = FALSE`, `useClimateForecast = FALSE`) over the large template, 1985 spin-up,
1986-01 to 2025-12.

```
Rscript backcast/Workflow_backcast_components.R [temperature] [ndvi] [flood] [inputs] [calibration] [forecast]
```

| Step | Output |
|---|---|
| `temperature` | `data/environmental_data/CDS/kriged_temperature_large`: ERA5-Land kriged to the large template, 1984-11 to 2025-12 (the earlier `kriged_temperature_v2` covers the park only). Needs CDS credentials and KrigR; slow |
| `ndvi` | Landsat 5/7 products rebuilt from the 1984+ downloads (2004+ versions kept in `<sensor>_2004`); `NDVI_large/harmonised_1984/{mean,sd,filled}.tif` |
| `flood` | `data/backcast/LAST_monthly_floodmaps_large` (monthly LAST maps from the WCS scenes) and `flood_monthly_500m_large_1984.tif` |
| `inputs` | `data/backcast/components/inputs_large`, `inputs_park`: one folder per component input |
| `calibration` | `data/backcast/rabbit_backcast_1985/calibration` (posterior, in-sample simulation) |
| `forecast` | `data/backcast/rabbit_backcast_1985/backcast` (component outputs, one density GeoTIFF per run), `rabbit_density_mean.tif`, `total_abundance.csv`, `agreement_with_census.csv` |

Each step skips what is on disk. `RABBIT_DOCKER=1` runs the components in their Docker images.
The products of the 2005-2025 lynx hindcast are left as they are.

## Files

- `Workflow_backcast_data.R`: climate before 2005 and of the current year, LAST flood maps from
  1984 and Landsat NDVI (phase 2), plus a coverage report. Every step skips what is on disk.
- `harmonise_landsat_to_s2.R`: linear adjustment of each Landsat sensor to Sentinel-2 (Landsat 5
  via Landsat 7), fitted on observed values only, with hold-out error in
  `data/backcast/harmonisation_fits.csv`; the combined series is gap-filled like the Sentinel-2 one.
- `Workflow_backcast_rabbit.R`: runs the model with posterior draws; writes maps, transect KAI,
  `total_abundance.csv/.png`, `backcast_vs_census.png` and `agreement_with_census.csv`.
- `functions/`: precipitation and flood download scripts, monthly flood composites, helpers.
  Temperature (`download_and_krig_temp()`) and NDVI (`build_monthly_ndvi()`) use the shared
  functions in `environmental_data/functions`.
- `ndvi_pipeline/`: copy of `NDVI-Pipeline-Donana-main_v2` with `--source landsat` and
  `--platforms`; Sentinel-2 behaviour is the same as the original.
