# Extending the large-template rabbit hindcast to 1990 / 1980: inputs and Earth Engine exports

*8 November 2026. Research note, nothing has been run on Earth Engine yet.*

The large-template hindcast (`backcast/rabbit_hindcast_large.R`) now runs 2005–2025 on
`data/environmental_data/template_raster_500_large.tif`. That grid is 189 × 173 cells at 500 m in
EPSG:25829, with extent x 669894.57–764394.57 and y 4069831.61–4156331.61. In WGS 84 that is
−7.0967 to −6.0085 E and 36.7374 to 37.5383 N, with 23,743 land cells. This note lists every
environmental input the model reads, where each one can come from back to 1980, and the Earth
Engine (GEE) exports needed. It also covers the local processing to do afterwards.

**Short answer.** A start in **1985** (1985 as spin‑up, results from 1986) is well supported. It
covers 1990, and it gives a few years before RHD arrived. Every input except NDVI before March 1984
and flooding before March/September 1984 is available. A start in **1980** is possible only with
much weaker NDVI (Landsat MSS or 5 km AVHRR) and no observed flooding for 1980–83. I do not
recommend it, except as a sensitivity run (section e).

The scripts are in `backcast/gee_1980/`:

| File | What it exports | Tasks |
|---|---|---|
| `00_scene_inventory.js` | CSV of every Landsat MSS/TM/ETM+ scene over the template, 1980–2003, with its clear share (**run first**) | 1 |
| `01_landsat_tm_ndvi_500m.js` | monthly 500 m NDVI mean/SD/valid fraction, Landsat 5 (1984–2003), Landsat 7 SLC‑on (1999–2003) and optionally Landsat 4 (1982–1993) | 20 + 5 (+12) |
| `02_jrc_water_500m.js` | monthly 500 m water fraction from JRC GSW, 1984‑03 to 2003 | 20 |
| `03_era5land_climate.js` | ERA5‑Land monthly t2m and precipitation, 1979–2025, native grid (optionally TerraClimate) | 2 (+2) |
| `04_avhrr_gimms_ndvi.js` | AVHRR CDR and GIMMS 3g monthly NDVI, 1981–2013 (only for a start before 1985) | 2 |
| `05_landsat_mss_ndvi_500m.js` | monthly 500 m NDVI from Landsat 2–5 MSS, 1980–1992 (only for a start before 1985) | 13 |
| `06_corine_landcover.js` | CORINE 1990 and 2000 at 100 m on the template grid | 2 |
| `split_yearly_gee_exports.R` | splits the yearly multi‑band exports into the per‑month files the R code reads (tested on a synthetic file) | – |

**Task count:** 48 tasks for a 1985/1990 start: 00, 01 (L5 and L7), 02, 03 and 06, of which 03 is
optional because the CDS route already works. A 1980 start needs about 29 more: L4, MSS, AVHRR/GIMMS
and TerraClimate.

---

## (a) What the model needs

From `rabbit_hindcast_large.R`, `rabbit/functions/Format_climate_data_for_rabbit.R`
(`calculate_BM_cDM`), `rabbit/functions/create_K_map_gaussian.R` and the
`environmental_data/create_large_*.R` scripts.

| Input | How the model uses it | Current source (large template) | Resolution / units | Current coverage | Needed from (start year Y) |
|---|---|---|---|---|---|
| **Temperature** | Monthly breeding probability (Tablado logistic in T, day length and its change, wet/dry), and the dry‑month test P < 2T | ERA5‑Land monthly mean 2 m temperature from CDS, kriged with KrigR on the GLO‑90 DEM to 0.0025° (`CDS/kriged_temperature_v2`), bilinear to 500 m. The western 2 km edge is filled by `fill_edges()` | monthly, K (converted to °C) | 2004‑11 → 2026 | **Nov Y−1** (the code stops otherwise) |
| **Precipitation** | Dry months (P < 2T) and the count of consecutive dry months | ERA5‑Land monthly means from CDS, `tp` (**m/day**, × days × 1000 gives mm/month), 0.1°, bilinear to 500 m | monthly | 2004‑01 → 2026‑07 | **Nov Y−1** |
| **NDVI mean** (500 m mean of 10 m NDVI) | Gaussian term of K (N_opt, N_sig) | Monthly max composite, aggregated on GEE (`ndvi_monthly_500m_gee.py`). Sentinel‑2 from 2015‑07; before that Landsat 8/7/5 adjusted to S2 with linear fits (`harmonise_landsat_to_s2.R`, `harmonise_area = "large"`). Gaps filled linearly per cell | monthly, unitless | 2004‑01 → 2026‑10 (`NDVI_large/harmonised/mean.tif`) | Nov Y−1 (gaps filled by `check_time_span_list(fill_gaps = TRUE)`) |
| **NDVI SD** (within‑cell SD of the 10 m pixels) | Gaussian term of K (S_opt, S_sig) | as NDVI mean (`.../harmonised/sd.tif`) | monthly | 2004‑01 → 2026‑10 | Nov Y−1 |
| **Flooding** (0/1) | K = 0 in flooded cells | In the park: LAST monthly flood maps (`LAST_monthly_floodmaps`, from the ICTS WCS). Elsewhere: JRC GSW (to 2016, water_fraction > 0) and Sentinel‑2 SCL (from 2017, > 0.05). No‑rabbit cells are always 1 (`create_large_flood_layer.R`) | monthly | 2004‑01 → 2025‑12 (LAST files locally from 2000; JRC exports from 2004) | Nov Y−1 |
| **Vegetation height** | Gaussian term of K (H_opt, H_sig) | **Static**: LiDAR in the park, ETH 2020 rescaled to the LiDAR elsewhere (`create_large_height_layer.R`) | m | single map (~2015–2020) | static (no historical product exists, see caveats) |
| **No‑rabbit mask** (towns, salines, water) | Forced to flooded (K = 0) every month | **Static**, CORINE 2018 (`create_large_landcover_layers.R`) | 500 m majority | 2018 | static; CLC 1990 is the better choice for the 1980s–90s |
| Day length | Breeding probability | Computed from latitude in `calculate_BM_cDM` | – | – | – |
| Template, park mask | Grid, park/outside split | Static | – | – | – |
| RHD arrival (not environmental) | Disease process | `rhd_settings`: RHD1 arrival ~1989 (`RHD1_arrival_year`), RHD2 2011‑10 | – | – | must be set for Y < 1989 |

The code does **not** read soil, land use other than the no‑rabbit mask, roads, or distance layers
(those are lynx inputs).

## (b) Recommended source per input, back to 1980

Dataset IDs, bands and date ranges were checked in the GEE data catalogue on 8 Nov 2026.

| Input | Recommended | GEE ID / bands | Coverage | Consistency with the current inputs |
|---|---|---|---|---|
| Temperature | **ERA5‑Land monthly, as now** (CDS + KrigR; the code supports any `start_date`). GEE as the alternative or cross‑check | `ECMWF/ERA5_LAND/MONTHLY_AGGR`, `temperature_2m` (K), 11 km | 1950‑02 → 2026‑10 | Same reanalysis as the calibration, so no bias correction is needed. 1979 onwards is the satellite‑assimilation era of ERA5. Check: pre‑2005 kriged vs raw ERA5 (as done for the relabelling) |
| Precipitation | **ERA5‑Land monthly, as now** (`download_precipitation_backcast.py --start-year 1979`) | same, `total_precipitation_sum` (**m per month**) | 1950 → | Same product. **Unit trap:** the GEE band is the monthly *total*, the CDS file is the mean *daily* total; divide the GEE values by days in the month. Negative or extreme GRIB‑packing values are possible: clamp at 0 |
| Independent climate check (optional) | TerraClimate | `IDAHO_EPSCOR/TERRACLIMATE`, `tmmx`, `tmmn` (× 0.1 °C), `pr` (mm), 4.6 km | 1958 → 2024 | Only to check ERA5‑Land trends in the 1980s; not a model input |
| NDVI 1984‑03 → 2003 | **Landsat 5 TM C2 L2**, same method as the existing 2004–2011 Landsat 5 files | `LANDSAT/LT05/C02/T1_L2`, `SR_B3` red, `SR_B4` NIR, `QA_PIXEL`; SR = DN × 0.0000275 − 0.2 | 1984‑03‑16 → 2012‑05 | Adjust to S2 through Landsat 7 as now. **New:** also export Landsat 7 1999‑06 → 2003‑05 (SLC‑on), so the L5→L7 fit gets 4 years of stripe‑free overlap (on the large template it now uses only 2004–2011, all SLC‑off) |
| NDVI, Landsat 4 TM (optional) | fills 1982–1993 gaps | `LANDSAT/LT04/C02/T1_L2`, same bands | 1982‑08 → 1993‑06 | Few scenes over Europe. Fit L4 → adjusted L5 (overlap 1984–1993) |
| NDVI 1980 → 1984‑02 (only for a 1980 start) | **Landsat MSS** TOA NDVI | `LANDSAT/LM02/C02/T1` (to 1982‑02), `LM03` (1978‑06 → 1983‑02): red `B5`, NIR `B7`; `LM04` (1982‑08 → 1992‑08), `LM05` (1984‑04 → 2013): red `B2`, NIR `B4`. TOA from `REFLECTANCE_MULT/ADD_BAND_x`; masks `QA_PIXEL` bits 0 and 3, `QA_RADSAT` | 1975 → 1992 | No SR, 60 m pixels, cloud bit only. Mean and SD are not comparable with TM. Fit MSS → adjusted L5 on the near‑simultaneous L4/L5 MSS + TM pairs (1982/84–1992). The SD term is the weakest part |
| NDVI, coarse support (only for a 1980 start, or to guide gap filling) | **NOAA CDR AVHRR NDVI v5** (GIMMS 3g as the alternative) | `NOAA/CDR/AVHRR/NDVI/V5`, `NDVI` × 0.0001, `QA` bits 1 cloud, 2 shadow, 3 water, 0.05°, daily; `NASA/GIMMS/3GV0` `ndvi`, `qa` ≤ 2, 1/12°, 15‑day | 1981‑06‑24 → 2013 (both) | 5.6/9 km. Use only as monthly **anomalies** applied to a per‑cell Landsat climatology, calibrated 1984–2011. Gives nothing for Jan–May 1981 or for 1980 |
| Flooding, park, 1984‑09 → | **LAST flood maps** (as now; not on GEE) | `backcast/functions/download_flooding_wcs_backcast.py` (ICTS WCS `inundacionv1`) → `make_monthly_floodmaps(min_year = 1984)` | 1984‑09‑10 → | Same source as the calibration |
| Flooding, outside the park and park months without LAST | **JRC GSW monthly history**, as now | `JRC/GSW1_4/MonthlyHistory`, `water` (0 no data, 1 land, 2 water), 30 m | 1984‑03 → 2021‑12 | Same product and threshold as 2004–2016. Before ~2000 it is Landsat 5 only, so more months are cloud‑ or data‑limited: check `valid_fraction` |
| Flooding before 1984‑03 | **no observed source** | (Landsat MSS water detection is not reliable over the marisma) | – | Use a flood climatology (section d) |
| No‑rabbit mask | **CORINE 1990** (Spain ~1987 imagery), CLC 2000 for the 1990s | `COPERNICUS/CORINE/V20/100m/1990`, `/2000`, band `landcover` (3‑digit codes) | 1990 (1986–1998 imagery), 2000 | GEE uses 3‑digit codes, the local CLC 2018 uses the 1–44 index, so recode first |
| Vegetation height | keep static | – | – | caveat |

## (c) Export scripts

### Common set‑up

- **Grid.** All 500 m exports use exactly the grid of `ndvi_monthly_500m_gee.py`: EPSG:32629
  (WGS 84 / UTM 29N, < 1 m from EPSG:25829), with the template origin rounded to whole 10 m:
  `x0 = 669890`, `y0 = 4156330`, 189 × 173 cells of 500 m, `crsTransform = [500,0,669890,0,-500,4156330]`.
  `build_monthly_ndvi_500m.R::on_template()` already shifts such files back onto the template
  (the shift is < 10 m). For EPSG:25829 output instead, set `OUTPUT_CRS = 'EPSG:25829'` and keep
  the transform; `on_template()` accepts either.
- **Method.** The NDVI and water exports are line‑by‑line ports of the existing Python scripts:
  monthly maximum NDVI on a 10 m default projection, then `reduceResolution` to mean, sample SD and
  valid fraction over the 50 × 50 sub‑pixels; scene cloud cover < 60 %; and QA_PIXEL bits 0, 1, 3,
  4 and 5 masked. New 1984–2003 files are therefore directly comparable with the existing 2004–2025
  files.
- **One task per year (Code Editor).** Every Code Editor export has to be started by hand in the
  Tasks tab. Monthly tasks would mean about 240 clicks for Landsat 5 alone, so the JS scripts export
  **one 36‑band GeoTIFF per year**, with bands `<YYYY>_<MM>_mean`, `_sd` and `_valid_fraction`.
  Months without scenes come out as −9999 with valid fraction 0. `split_yearly_gee_exports.R` turns
  them into the per‑month files of the Python workflow.
- **Python alternative, no clicking.** The existing scripts already take any start date:

  ```
  python backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py --project <ee-project> \
     --template data/environmental_data/template_raster_500_large.tif \
     --sensor landsat-5 --start 1984-03 --end 2003-12 --drive-folder ndvi500_large_landsat-5
  python backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py ... --sensor landsat-7 --start 1999-06 --end 2003-05 \
     --drive-folder ndvi500_large_landsat-7
  python backcast/ndvi_pipeline/height_water_gee.py water --project <ee-project> \
     --template data/environmental_data/template_raster_500_large.tif \
     --source jrc --start 1984-03 --end 2003-12 --drive-folder water500_large_jrc
  ```

  These submit about 240, 48 and 238 tasks. GEE queues them and runs a few at a time. The Python
  script cannot do Landsat 4 or MSS; the JS scripts (01 with `landsat-4`, and 05) cover those.
- **Limits.** If a yearly NDVI task fails with "User memory limit exceeded" or "Computation timed
  out", lower `shardSize` from 40 to 20, or split the year into two 6‑month images. A 36‑band
  shard is 12 times the work of one Python task.
- **Landsat 7 SLC‑off.** From 31 May 2003 every Landsat 7 scene has about 22 % striped no‑data.
  Script 01 keeps Landsat 7 to 1999‑05 → 2003‑05. The existing 2004–2022 Landsat 7 series
  (SLC‑off) stays as it is: in the harmonised series Landsat 7 is used last, after L8 and L5.

The full scripts follow. They are identical to the files in `backcast/gee_1980/`.


### `00_scene_inventory.js`

```javascript
// =====================================================================================
// 00_scene_inventory.js  -  Landsat scene inventory over the large Doñana template, 1980-2003
//
// Run this FIRST. It lists every Landsat MSS/TM/ETM+ scene that touches the large (lynx)
// template, with its scene cloud cover and the share of the template that is clear (QA_PIXEL),
// so you can see which months of the 1980s/1990s have usable NDVI before submitting the
// heavy NDVI exports.
//
// Output: one CSV in Google Drive (folder gee_backcast_1980), one row per scene:
//   collection, sensor, date, year, month, CLOUD_COVER, clear_frac (share of the template
//   bounding box with a clear, non-fill pixel), WRS_PATH, WRS_ROW
// One export task. Also prints a chart of TM scenes per month in the Console.
// =====================================================================================

var OUTPUT_CRS = 'EPSG:32629';   // as backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py
// Large template (data/environmental_data/template_raster_500_large.tif, 189 x 173 cells, 500 m),
// origin rounded to whole 10 m as in ndvi_monthly_500m_gee.py
var X0 = 669890, Y0 = 4156330, RES = 500, NCOL = 189, NROW = 173;
var region = ee.Geometry.Rectangle([X0, Y0 - NROW * RES, X0 + NCOL * RES, Y0],
                                   ee.Projection(OUTPUT_CRS), false);
var START = '1980-01-01', END = '2004-01-01';
var DRIVE_FOLDER = 'gee_backcast_1980';

// collection id, sensor label, QA bits that make a pixel unusable
var TM_BAD  = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4) | (1 << 5);  // fill, dilated cloud, cloud, shadow, snow
var MSS_BAD = (1 << 0) | (1 << 3);                                   // fill, cloud (MSS has no shadow bit)
var COLLECTIONS = [
  ['LANDSAT/LM02/C02/T1',    'mss-2', MSS_BAD],
  ['LANDSAT/LM03/C02/T1',    'mss-3', MSS_BAD],
  ['LANDSAT/LM04/C02/T1',    'mss-4', MSS_BAD],
  ['LANDSAT/LM05/C02/T1',    'mss-5', MSS_BAD],
  ['LANDSAT/LT04/C02/T1_L2', 'landsat-4', TM_BAD],
  ['LANDSAT/LT05/C02/T1_L2', 'landsat-5', TM_BAD],
  ['LANDSAT/LT05/C02/T2_L2', 'landsat-5-T2', TM_BAD],   // Tier 2: poorer geolocation, may fill 1980s gaps
  ['LANDSAT/LE07/C02/T1_L2', 'landsat-7', TM_BAD]
];

function inventory(id, sensor, bad) {
  return ee.ImageCollection(id).filterBounds(region).filterDate(START, END).map(function (img) {
    var clear = img.select('QA_PIXEL').bitwiseAnd(bad).eq(0).unmask(0);
    var frac = clear.reduceRegion({reducer: ee.Reducer.mean(), geometry: region, scale: 300,
                                   crs: OUTPUT_CRS, maxPixels: 1e8, bestEffort: true}).values().get(0);
    var d = ee.Date(img.get('system:time_start'));
    return ee.Feature(null, {
      collection: id, sensor: sensor, date: d.format('YYYY-MM-dd'),
      year: d.get('year'), month: d.get('month'),
      CLOUD_COVER: img.get('CLOUD_COVER'), clear_frac: frac,
      WRS_PATH: img.get('WRS_PATH'), WRS_ROW: img.get('WRS_ROW')
    });
  });
}

var all = ee.FeatureCollection(COLLECTIONS.map(function (c) { return inventory(c[0], c[1], c[2]); })).flatten();

Export.table.toDrive({
  collection: all, description: 'landsat_scene_inventory_1980_2003', folder: DRIVE_FOLDER,
  fileNamePrefix: 'landsat_scene_inventory_1980_2003', fileFormat: 'CSV',
  selectors: ['collection', 'sensor', 'date', 'year', 'month', 'CLOUD_COVER', 'clear_frac', 'WRS_PATH', 'WRS_ROW']
});

// Quick look: Landsat 5 TM T1 scenes per year in the Console
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(region).filterDate('1984-01-01', '2004-01-01');
var years = ee.List.sequence(1984, 2003);
print('Landsat 5 TM T1 scenes per year (all cloud covers)',
      ee.Dictionary.fromLists(years.map(function (y) { return ee.Number(y).format('%d'); }),
                              years.map(function (y) { return l5.filter(ee.Filter.calendarRange(y, y, 'year')).size(); })));
Map.centerObject(region, 9);
Map.addLayer(region, {color: 'red'}, 'large template');
```

### `01_landsat_tm_ndvi_500m.js`

```javascript
// =====================================================================================
// 01_landsat_tm_ndvi_500m.js  -  monthly 500 m NDVI mean / SD / valid fraction from
//                                Landsat 4/5 TM and Landsat 7 ETM+ (Collection 2 Level-2)
//
// Code Editor port of backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py (same method, so the
// output joins the existing 2004-2025 files in data/environmental_data/ndvi500_large/):
//   1. all scenes of the month with scene CLOUD_COVER < 60, masked per pixel with QA_PIXEL
//      bits 0 fill, 1 dilated cloud, 3 cloud, 4 cloud shadow, 5 snow;
//      SR = DN * 0.0000275 - 0.2; NDVI = (NIR - red) / (NIR + red), clamped to [-1, 1]
//   2. monthly composite = maximum NDVI per pixel, on a 10 m grid (30 m pixels kept, nearest)
//   3. per 500 m template cell: mean, sample SD (n - 1) and valid fraction (share of the
//      50 x 50 10 m pixels with a value)
//
// Difference from the Python script: to keep the number of tasks low in the Code Editor
// (each task has to be started by hand), ONE task per year exports a 36-band image:
//   ndvi500_<sensor>_<YYYY>.tif   bands <YYYY>_<MM>_mean, <YYYY>_<MM>_sd, <YYYY>_<MM>_valid_fraction
//   -9999 = no data; months without any scene have mean/sd -9999 and valid_fraction 0.
// split_yearly_gee_exports.R (see gee_backcast_1980_requirements.md) turns them into the
// per-month ndvi500_<sensor>_<YYYY-MM>.tif files that build_monthly_ndvi_500m() reads.
//
// Alternative with one task per month and no clicking: the Python script itself, e.g.
//   python backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py --project <ee-project> \
//     --template data/environmental_data/template_raster_500_large.tif \
//     --sensor landsat-5 --start 1984-03 --end 2003-12 --drive-folder ndvi500_large_landsat-5
//
// Landsat 7: SLC failed on 31 May 2003; afterwards every scene has ~22% striped no-data.
// Export Landsat 7 1999-06 .. 2003-05 here: SLC-on overlap with Landsat 5 for the harmonisation.
// Landsat 4 TM: 1982-08 .. 1993-06 (very few scenes over Europe), only useful for the 1980s.
//
// Usage: set SENSOR, FIRST_YEAR, LAST_YEAR, press Run, then start each task in the Tasks tab.
// Suggested runs:
//   landsat-5  1984-2003   (20 tasks)
//   landsat-7  1999-2003   (5 tasks; months after May 2003 are skipped)
//   landsat-4  1982-1993   (12 tasks, optional, for 1980s gaps)
// =====================================================================================

var SENSOR = 'landsat-5';
var FIRST_YEAR = 1984, LAST_YEAR = 2003;
var CLOUD_COVER = 60;   // scene cloud cover limit, %, as the Python default
var DRIVE_FOLDER = 'ndvi500_large_' + SENSOR + '_1984';

var SENSORS = {
  'landsat-4': {collection: 'LANDSAT/LT04/C02/T1_L2', red: 'SR_B3', nir: 'SR_B4', first: '1982-08-01', last: '1993-07-01'},
  'landsat-5': {collection: 'LANDSAT/LT05/C02/T1_L2', red: 'SR_B3', nir: 'SR_B4', first: '1984-03-01', last: '2012-06-01'},
  'landsat-7': {collection: 'LANDSAT/LE07/C02/T1_L2', red: 'SR_B3', nir: 'SR_B4', first: '1999-05-01', last: '2003-06-01'}  // SLC-on only
};
var spec = SENSORS[SENSOR];

// ---- template grid (template_raster_500_large.tif, origin rounded to 10 m) -------------
var OUTPUT_CRS = 'EPSG:32629', FINE = 10, NODATA = -9999;
var X0 = 669890, Y0 = 4156330, RES = 500, NCOL = 189, NROW = 173, FACT = RES / FINE;
var region = ee.Geometry.Rectangle([X0, Y0 - NROW * RES, X0 + NCOL * RES, Y0],
                                   ee.Projection(OUTPUT_CRS), false);
var FINE_T   = [FINE, 0, X0, 0, -FINE, Y0];
var COARSE_T = [RES, 0, X0, 0, -RES, Y0];

var BAD = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4) | (1 << 5);
var SCALE = 0.0000275, OFFSET = -0.2;

function ndvi(img) {
  var valid = img.select('QA_PIXEL').bitwiseAnd(BAD).eq(0);
  var red = img.select(spec.red).multiply(SCALE).add(OFFSET);
  var nir = img.select(spec.nir).multiply(SCALE).add(OFFSET);
  var den = nir.add(red);
  return nir.subtract(red).divide(den).updateMask(den.neq(0)).clamp(-1, 1)
            .updateMask(valid).rename('ndvi');
}

var EMPTY = ee.Image.constant([NODATA, NODATA, 0]).rename(['mean', 'sd', 'valid_fraction']).toFloat();

function monthStats(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var col = ee.ImageCollection(spec.collection).filterBounds(region)
              .filterDate(start, start.advance(1, 'month'))
              .filterDate(spec.first, spec.last)
              .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_COVER));
  var reducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true)
                                 .combine(ee.Reducer.count(), null, true);
  var comp = col.map(ndvi).max().setDefaultProjection({crs: OUTPUT_CRS, crsTransform: FINE_T});
  var st = comp.reduceResolution({reducer: reducer, maxPixels: FACT * FACT + 100})
               .reproject({crs: OUTPUT_CRS, crsTransform: COARSE_T});
  var n = st.select('ndvi_count');
  // GEE stdDev is the population SD; R's sd() is the sample SD
  var sd = st.select('ndvi_stdDev').multiply(n.divide(n.subtract(1)).sqrt()).updateMask(n.gt(1));
  var out = ee.Image.cat([st.select('ndvi_mean').rename('mean'), sd.rename('sd'),
                          n.divide(FACT * FACT).unmask(0).rename('valid_fraction')])
              .toFloat().unmask(NODATA);
  var img = ee.Image(ee.Algorithms.If(col.size().gt(0), out, EMPTY));
  var p = year + '_' + (month < 10 ? '0' : '') + month + '_';
  return img.rename([p + 'mean', p + 'sd', p + 'valid_fraction']);
}

for (var y = FIRST_YEAR; y <= LAST_YEAR; y++) {
  var months = [];
  for (var m = 1; m <= 12; m++) months.push(monthStats(y, m));
  var name = 'ndvi500_' + SENSOR + '_' + y;
  Export.image.toDrive({
    image: ee.Image.cat(months), description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
    region: region, crs: OUTPUT_CRS, crsTransform: COARSE_T, maxPixels: 1e9, fileFormat: 'GeoTIFF',
    // shards of 40 x 40 output cells = 2000 x 2000 pixels at 10 m (the default 256 fails with
    // "Reprojection output too large")
    shardSize: 40
  });
}
print('Submitted (start them in the Tasks tab):', SENSOR, FIRST_YEAR + '-' + LAST_YEAR);
```

### `02_jrc_water_500m.js`

```javascript
// =====================================================================================
// 02_jrc_water_500m.js  -  monthly surface water fraction per 500 m cell from the JRC Global
//                          Surface Water monthly history (Landsat 5/7/8, 30 m), 1984-03 .. 2003-12
//
// Code Editor port of `height_water_gee.py water --source jrc` (same method, so the output joins
// data/environmental_data/water_large/jrc/water500_jrc_<YYYY-MM>.tif, which start in 2004-01):
//   per 30 m pixel: 1 = water, 0 = observed not water, masked = no observation (JRC value 0);
//   per 500 m cell (via a 10 m grid, nearest): water_fraction = mean, valid_fraction = share of
//   the 50 x 50 10 m pixels observed.
// create_large_flood_layer.R uses water_fraction > 0 = flooded, outside the park and in park
// months without a LAST flood map.
//
// One task per year: water500_jrc_<YYYY>.tif with bands <YYYY>_<MM>_water_fraction and
// <YYYY>_<MM>_valid_fraction (-9999 = no data; months before 1984-03 are all -9999 / 0).
// Split into monthly files with split_yearly_gee_exports.R.
// Python alternative (one task per month):
//   python backcast/ndvi_pipeline/height_water_gee.py water --project <ee-project> \
//     --template data/environmental_data/template_raster_500_large.tif \
//     --source jrc --start 1984-03 --end 2003-12 --drive-folder water500_large_jrc
// =====================================================================================

var FIRST_YEAR = 1984, LAST_YEAR = 2003;
var DRIVE_FOLDER = 'water500_large_jrc_1984';

var OUTPUT_CRS = 'EPSG:32629', FINE = 10, NODATA = -9999;
var X0 = 669890, Y0 = 4156330, RES = 500, NCOL = 189, NROW = 173, FACT = RES / FINE;
var region = ee.Geometry.Rectangle([X0, Y0 - NROW * RES, X0 + NCOL * RES, Y0],
                                   ee.Projection(OUTPUT_CRS), false);
var FINE_T   = [FINE, 0, X0, 0, -FINE, Y0];
var COARSE_T = [RES, 0, X0, 0, -RES, Y0];

var JRC = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory');
var EMPTY = ee.Image.constant([NODATA, 0]).rename(['water_fraction', 'valid_fraction']).toFloat();

function monthWater(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var col = JRC.filterDate(start, start.advance(1, 'month'));
  var w = ee.Image(col.first()).select('water');
  var water = w.eq(2).updateMask(w.gt(0)).rename('water').toFloat()
               .setDefaultProjection({crs: OUTPUT_CRS, crsTransform: FINE_T});
  var st = water.reduceResolution({
      reducer: ee.Reducer.mean().combine(ee.Reducer.count(), null, true),
      maxPixels: FACT * FACT + 100})
    .reproject({crs: OUTPUT_CRS, crsTransform: COARSE_T});
  var out = ee.Image.cat([st.select('water_mean').rename('water_fraction'),
                          st.select('water_count').divide(FACT * FACT).unmask(0).rename('valid_fraction')])
              .toFloat().unmask(NODATA);
  var img = ee.Image(ee.Algorithms.If(col.size().gt(0), out, EMPTY));
  var p = year + '_' + (month < 10 ? '0' : '') + month + '_';
  return img.rename([p + 'water_fraction', p + 'valid_fraction']);
}

for (var y = FIRST_YEAR; y <= LAST_YEAR; y++) {
  var months = [];
  for (var m = 1; m <= 12; m++) months.push(monthWater(y, m));
  var name = 'water500_jrc_' + y;
  Export.image.toDrive({
    image: ee.Image.cat(months), description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
    region: region, crs: OUTPUT_CRS, crsTransform: COARSE_T, maxPixels: 1e9, fileFormat: 'GeoTIFF',
    shardSize: 40
  });
}
print('Submitted (start them in the Tasks tab): JRC water', FIRST_YEAR + '-' + LAST_YEAR);
```

### `03_era5land_climate.js`

```javascript
// =====================================================================================
// 03_era5land_climate.js  -  monthly ERA5-Land 2 m temperature and precipitation, 1950 .. now,
//                            on the native ERA5-Land grid (0.1 deg), around the large template
//
// The rabbit model uses
//   temperature: ERA5-Land monthly mean 2 m temperature (K), kriged to ~250 m with the DEM as
//                covariate (environmental_data/functions/Download_and_Krig_temperature.R)
//   precipitation: ERA5-Land monthly precipitation, bilinear to the 500 m template
// Both currently come from the Copernicus CDS ("reanalysis-era5-land-monthly-means"), which
// already goes back to 1950: the CDS route (download_and_krig_temp(start_date = ...) and
// backcast/functions/download_precipitation_backcast.py --start-year 1979) is the preferred,
// identical-to-calibration way. This script is the Earth Engine alternative / cross-check.
//
// !! UNITS: GEE's ECMWF/ERA5_LAND/MONTHLY_AGGR band total_precipitation_sum is the MONTHLY TOTAL
//    in m. The CDS monthly-means file the model reads (tp) is the MEAN DAILY total in m/day, and
//    calculate_BM_cDM() multiplies it by days_in_month * 1000. Divide the GEE values by the
//    number of days in the month before feeding them to the model (or the model would see
//    ~30x the rain). temperature_2m is the monthly mean in K in both.
//
// Output (Google Drive folder gee_backcast_1980), one task per variable:
//   era5land_t2m_monthly_<FIRST>_<LAST>.tif   bands YYYY_MM, K
//   era5land_tp_monthly_<FIRST>_<LAST>.tif    bands YYYY_MM, m per month
// Region: template bounding box + 0.5 deg on each side (as the kriging download area).
//
// Optional (EXPORT_TERRACLIMATE = true): TerraClimate (4 km, 1958 .. 2024), monthly mean of
// tmax/tmin (deg C) and precipitation (mm/month) - only as an independent check of ERA5-Land
// trends in the 1980s, not as model input.
// =====================================================================================

var FIRST = '1979-01-01', LAST = '2026-01-01';   // end date exclusive
var DRIVE_FOLDER = 'gee_backcast_1980';
var EXPORT_TERRACLIMATE = false;

// large template in WGS 84: -7.09665 .. -6.00849 E, 36.73741 .. 37.53833 N, + 0.5 deg
var region = ee.Geometry.Rectangle([-7.6, 36.2, -5.5, 38.1], null, false);

var era = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').filterDate(FIRST, LAST);
var names = era.aggregate_array('system:time_start').map(function (t) {
  return ee.Date(t).format('YYYY_MM');
});
var proj = era.first().projection().getInfo();   // native 0.1 deg grid, no resampling
print('ERA5-Land native projection', proj, 'months', names.size());

function exportVar(band, short) {
  var img = era.select(band).toBands().rename(names).toFloat();
  var name = 'era5land_' + short + '_monthly_' + FIRST.slice(0, 4) + '_' + (Number(LAST.slice(0, 4)) - 1);
  Export.image.toDrive({image: img, description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
                        region: region, crs: proj.crs, crsTransform: proj.transform,
                        maxPixels: 1e9, fileFormat: 'GeoTIFF'});
}
exportVar('temperature_2m', 't2m');
exportVar('total_precipitation_sum', 'tp');

if (EXPORT_TERRACLIMATE) {
  var tc = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate(FIRST, LAST);
  var tcNames = tc.aggregate_array('system:time_start').map(function (t) { return ee.Date(t).format('YYYY_MM'); });
  var tcProj = tc.first().projection().getInfo();
  var tmean = tc.map(function (i) { return i.select('tmmx').add(i.select('tmmn')).multiply(0.05).rename('tmean'); })
                .toBands().rename(tcNames).toFloat();
  var pr = tc.select('pr').toBands().rename(tcNames).toFloat();
  [[tmean, 'terraclimate_tmean_degC'], [pr, 'terraclimate_pr_mm']].forEach(function (x) {
    Export.image.toDrive({image: x[0], description: x[1], folder: DRIVE_FOLDER, fileNamePrefix: x[1],
                          region: region, crs: tcProj.crs, crsTransform: tcProj.transform,
                          maxPixels: 1e9, fileFormat: 'GeoTIFF'});
  });
}
```

### `04_avhrr_gimms_ndvi.js`

```javascript
// =====================================================================================
// 04_avhrr_gimms_ndvi.js  -  coarse monthly NDVI (AVHRR) around the large template, 1981 .. 2013
//
// Only needed for a start before 1985 (or to guide gap filling in sparse Landsat 5 years).
// AVHRR pixels (5.6 km / 9 km) are 10-20 template cells wide: they carry the month-to-month
// and year-to-year signal, not the within-cell SD. Use them as anomalies on top of a Landsat
// per-cell climatology (see gee_backcast_1980_requirements.md, post-processing), calibrated in
// the Landsat 5 overlap (1984-2011).
//
//   NOAA/CDR/AVHRR/NDVI/V5   daily, 0.05 deg, 1981-06-24 .. 2013-12-31, NDVI * 0.0001,
//                            QA bit 1 cloud, bit 2 cloud shadow, bit 3 water
//   NASA/GIMMS/3GV0          15-day, 1/12 deg, 1981-07 .. 2013-12, ndvi (-1..1),
//                            qa 1-2 good, 3-6 interpolated/seasonal profile, 7 missing
//
// Monthly maximum-value composite of good pixels (as the Landsat/Sentinel-2 composites are
// monthly maxima), plus the number of good observations. Native grids, no resampling.
// Output (folder gee_backcast_1980), two tasks:
//   avhrr_cdr_ndvi_monthly_1981_2013.tif   bands YYYY_MM_ndvi, YYYY_MM_n  (-9999 = no data)
//   gimms3g_ndvi_monthly_1981_2013.tif     bands YYYY_MM_ndvi, YYYY_MM_n
// =====================================================================================

var DRIVE_FOLDER = 'gee_backcast_1980';
var NODATA = -9999;
var region = ee.Geometry.Rectangle([-7.2, 36.6, -5.9, 37.65], null, false);  // template + ~0.1 deg

// toBands() prefixes band names with the image index; build the names explicitly instead
function stack(col, prep, first, last) {
  var n = ee.Date(last).difference(ee.Date(first), 'month').round();
  var idx = ee.List.sequence(0, n.subtract(1));
  var imgs = idx.map(function (i) {
    var s = ee.Date(first).advance(i, 'month');
    var c = col.filterDate(s, s.advance(1, 'month')).map(prep);
    return ee.Image(ee.Algorithms.If(c.size().gt(0),
      ee.Image.cat([c.max().rename('ndvi'), c.count().rename('n')]),
      ee.Image.constant([NODATA, 0]).rename(['ndvi', 'n']))).toFloat().unmask(NODATA);
  });
  var names = idx.map(function (i) {
    var p = ee.Date(first).advance(i, 'month').format('YYYY_MM');
    return [ee.String(p).cat('_ndvi'), ee.String(p).cat('_n')];
  }).flatten();
  return ee.ImageCollection.fromImages(imgs).toBands().rename(names);
}

// --- NOAA CDR AVHRR NDVI V5 -----------------------------------------------------------
var cdr = ee.ImageCollection('NOAA/CDR/AVHRR/NDVI/V5').filterBounds(region);
var cdrBad = (1 << 1) | (1 << 2) | (1 << 3);
function cdrPrep(i) {
  return i.select('NDVI').multiply(0.0001)
          .updateMask(i.select('QA').bitwiseAnd(cdrBad).eq(0))
          .updateMask(i.select('NDVI').gt(-9998));
}
var cdrProj = cdr.first().projection().getInfo();
Export.image.toDrive({image: stack(cdr, cdrPrep, '1981-06-01', '2014-01-01'),
                      description: 'avhrr_cdr_ndvi_monthly_1981_2013', folder: DRIVE_FOLDER,
                      fileNamePrefix: 'avhrr_cdr_ndvi_monthly_1981_2013', region: region,
                      crs: cdrProj.crs, crsTransform: cdrProj.transform, maxPixels: 1e9});

// --- GIMMS NDVI3g -----------------------------------------------------------------------
var gimms = ee.ImageCollection('NASA/GIMMS/3GV0').filterBounds(region);
function gimmsPrep(i) {
  return i.select('ndvi').updateMask(i.select('qa').lte(2));
}
var gProj = gimms.first().projection().getInfo();
Export.image.toDrive({image: stack(gimms, gimmsPrep, '1981-07-01', '2014-01-01'),
                      description: 'gimms3g_ndvi_monthly_1981_2013', folder: DRIVE_FOLDER,
                      fileNamePrefix: 'gimms3g_ndvi_monthly_1981_2013', region: region,
                      crs: gProj.crs, crsTransform: gProj.transform, maxPixels: 1e9});
```

### `05_landsat_mss_ndvi_500m.js`

```javascript
// =====================================================================================
// 05_landsat_mss_ndvi_500m.js  -  monthly 500 m NDVI mean / SD / valid fraction from Landsat
//                                 MSS (Collection 2 Tier 1, Level-1), 1980 .. 1992
//
// Optional: the only 30-60 m NDVI before Landsat 5 TM (March 1984). Check the inventory
// (00_scene_inventory.js) first - MSS scenes over Doñana in 1980-1983 may be few.
//
// MSS has no surface-reflectance product. NDVI is computed from top-of-atmosphere reflectance
// (DN * REFLECTANCE_MULT_BAND_x + REFLECTANCE_ADD_BAND_x; the sun-elevation correction cancels
// in the ratio) with red = 0.6-0.7 um and NIR = 0.8-1.1 um (the closest to TM bands 3 and 4):
//   Landsat 1-3 MSS: red B5, NIR B7        Landsat 4-5 MSS: red B2, NIR B4
// Masks: QA_PIXEL bit 0 (fill) and bit 3 (cloud; there is no cloud-shadow flag),
// QA_RADSAT != 0 (saturated). Native resolution 60 m (delivered on a 30 m grid).
// TOA, 60 m, and a cruder cloud mask: the values are NOT comparable with TM. Harmonise them
// against Landsat 5 TM (and Landsat 4 TM) using the near-simultaneous MSS + TM acquisitions on
// the same platforms (Landsat 4: 1982-1992, Landsat 5: 1984-1992) - hence the export to 1992.
//
// Same aggregation as 01_landsat_tm_ndvi_500m.js; all MSS platforms are merged into one series
// ('mss'). One task per year: ndvi500_mss_<YYYY>.tif, bands <YYYY>_<MM>_{mean,sd,valid_fraction}.
// =====================================================================================

var FIRST_YEAR = 1980, LAST_YEAR = 1992;
var CLOUD_COVER = 60;
var DRIVE_FOLDER = 'ndvi500_large_mss';

var PLATFORMS = [
  {collection: 'LANDSAT/LM02/C02/T1', red: 'B5', nir: 'B7'},
  {collection: 'LANDSAT/LM03/C02/T1', red: 'B5', nir: 'B7'},
  {collection: 'LANDSAT/LM04/C02/T1', red: 'B2', nir: 'B4'},
  {collection: 'LANDSAT/LM05/C02/T1', red: 'B2', nir: 'B4'}
];

var OUTPUT_CRS = 'EPSG:32629', FINE = 10, NODATA = -9999;
var X0 = 669890, Y0 = 4156330, RES = 500, NCOL = 189, NROW = 173, FACT = RES / FINE;
var region = ee.Geometry.Rectangle([X0, Y0 - NROW * RES, X0 + NCOL * RES, Y0],
                                   ee.Projection(OUTPUT_CRS), false);
var FINE_T   = [FINE, 0, X0, 0, -FINE, Y0];
var COARSE_T = [RES, 0, X0, 0, -RES, Y0];
var BAD = (1 << 0) | (1 << 3);

function mssNdvi(p) {
  var rn = p.red.slice(1), nn = p.nir.slice(1);
  return function (img) {
    var red = img.select(p.red).multiply(ee.Number(img.get('REFLECTANCE_MULT_BAND_' + rn)))
                 .add(ee.Number(img.get('REFLECTANCE_ADD_BAND_' + rn)));
    var nir = img.select(p.nir).multiply(ee.Number(img.get('REFLECTANCE_MULT_BAND_' + nn)))
                 .add(ee.Number(img.get('REFLECTANCE_ADD_BAND_' + nn)));
    var valid = img.select('QA_PIXEL').bitwiseAnd(BAD).eq(0).and(img.select('QA_RADSAT').eq(0));
    var den = nir.add(red);
    return nir.subtract(red).divide(den).updateMask(den.neq(0)).clamp(-1, 1)
              .updateMask(valid).rename('ndvi').toFloat();
  };
}

var mss = PLATFORMS.map(function (p) {
  return ee.ImageCollection(p.collection).filterBounds(region)
           .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_COVER)).map(mssNdvi(p));
}).reduce(function (a, b) { return a.merge(b); });

var EMPTY = ee.Image.constant([NODATA, NODATA, 0]).rename(['mean', 'sd', 'valid_fraction']).toFloat();

function monthStats(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var col = mss.filterDate(start, start.advance(1, 'month'));
  var reducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true)
                                 .combine(ee.Reducer.count(), null, true);
  var st = col.max().setDefaultProjection({crs: OUTPUT_CRS, crsTransform: FINE_T})
              .reduceResolution({reducer: reducer, maxPixels: FACT * FACT + 100})
              .reproject({crs: OUTPUT_CRS, crsTransform: COARSE_T});
  var n = st.select('ndvi_count');
  var sd = st.select('ndvi_stdDev').multiply(n.divide(n.subtract(1)).sqrt()).updateMask(n.gt(1));
  var out = ee.Image.cat([st.select('ndvi_mean').rename('mean'), sd.rename('sd'),
                          n.divide(FACT * FACT).unmask(0).rename('valid_fraction')])
              .toFloat().unmask(NODATA);
  var img = ee.Image(ee.Algorithms.If(col.size().gt(0), out, EMPTY));
  var pfx = year + '_' + (month < 10 ? '0' : '') + month + '_';
  return img.rename([pfx + 'mean', pfx + 'sd', pfx + 'valid_fraction']);
}

for (var y = FIRST_YEAR; y <= LAST_YEAR; y++) {
  var months = [];
  for (var m = 1; m <= 12; m++) months.push(monthStats(y, m));
  var name = 'ndvi500_mss_' + y;
  Export.image.toDrive({
    image: ee.Image.cat(months), description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
    region: region, crs: OUTPUT_CRS, crsTransform: COARSE_T, maxPixels: 1e9, fileFormat: 'GeoTIFF',
    shardSize: 40
  });
}
print('Submitted (start them in the Tasks tab): MSS', FIRST_YEAR + '-' + LAST_YEAR);
```

### `06_corine_landcover.js`

```javascript
// =====================================================================================
// 06_corine_landcover.js  -  CORINE Land Cover 1990 and 2000 (100 m) on a 100 m grid aligned with
//                            the large template, for a period-specific rabbit no-habitat mask
//
// create_large_landcover_layers.R makes rabbit_no_habitat_500_large.tif (towns, salines, water:
// K = 0) from CLC 2018 and holds it fixed. For a 1980s/1990s backcast CLC 1990 (Spain: imagery
// of ~1987) is closer: greenhouses/berry farms, Matalascañas and other built-up areas grew a lot
// after 1990. The local CLC 2018 file uses the EEA raster index (1-44, 48); GEE uses the 3-digit
// CLC codes (111 .. 523) - convert in R before reusing the class lists of
// create_large_landcover_layers.R (index i <-> code: 1=111, 2=112, 3=121, ..., 44=523).
// The 1990 and 2000 rasters can also be downloaded directly from land.copernicus.eu.
//
// Output (folder gee_backcast_1980): corine<YEAR>_100m_large.tif, band landcover (3-digit code),
// on the template grid at 100 m (EPSG:32629, origin rounded to 10 m as the other exports).
// Two tasks.
// =====================================================================================

var DRIVE_FOLDER = 'gee_backcast_1980';
var OUTPUT_CRS = 'EPSG:32629';
var X0 = 669890, Y0 = 4156330, RES = 500, NCOL = 189, NROW = 173;
var region = ee.Geometry.Rectangle([X0, Y0 - NROW * RES, X0 + NCOL * RES, Y0],
                                   ee.Projection(OUTPUT_CRS), false);

['1990', '2000'].forEach(function (yr) {
  var img = ee.Image('COPERNICUS/CORINE/V20/100m/' + yr).select('landcover').toUint16();
  var name = 'corine' + yr + '_100m_large';
  Export.image.toDrive({image: img, description: name, folder: DRIVE_FOLDER, fileNamePrefix: name,
                        region: region, crs: OUTPUT_CRS, crsTransform: [100, 0, X0, 0, -100, Y0],
                        maxPixels: 1e9, fileFormat: 'GeoTIFF'});
});
```

### `split_yearly_gee_exports.R`

```r
#
# Split the yearly multi-band exports of gee_1980/01, 02 and 05 (*.js) into the per-month files
# that the existing R code reads:
#   ndvi500_<sensor>_<YYYY>.tif  ->  ndvi500_<sensor>_<YYYY-MM>.tif  (bands mean, sd, valid_fraction;
#                                    read by build_monthly_ndvi_500m())
#   water500_jrc_<YYYY>.tif      ->  water500_jrc_<YYYY-MM>.tif      (bands water_fraction,
#                                    valid_fraction; read by create_large_flood_layer.R)
# Months without any observation (valid_fraction 0 everywhere) are not written, as the Python
# exporters do not export months without scenes (they become gaps, filled later).
#
#   Rscript backcast/gee_1980/split_yearly_gee_exports.R <folder with yearly files> <output folder>

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "^(ndvi500|water500)_.+_\\d{4}(-\\d+-\\d+)?\\.tif$", full.names = TRUE)
for (f in files) {
  r <- rast(f)
  # GEE splits large exports into tiles (name-0000000000-0000000000.tif); the template is small
  # enough for one file, but merge if needed before running this.
  stem <- sub("_(\\d{4})(-\\d+-\\d+)?\\.tif$", "", basename(f))
  ym <- unique(sub("^(\\d{4})_(\\d{2})_.*$", "\\1-\\2", names(r)))
  for (m in ym) {
    lyr <- r[[startsWith(names(r), paste0(sub("-", "_", m), "_"))]]
    names(lyr) <- sub("^\\d{4}_\\d{2}_", "", names(lyr))
    if (all(values(lyr[["valid_fraction"]], mat = FALSE) %in% c(0, NA))) next
    lyr <- lyr[[if (startsWith(stem, "ndvi500")) c("mean", "sd", "valid_fraction")
                else c("water_fraction", "valid_fraction")]]
    writeRaster(lyr, file.path(out_dir, sprintf("%s_%s.tif", stem, m)), overwrite = TRUE,
                NAflag = -9999)
  }
  cat(basename(f), "->", length(ym), "months checked\n")
}
```

---

## (d) Local post-processing after the downloads

Run everything from the project root. The step order follows the existing workflow, and every
script listed here exists already; most only need their hard‑coded dates or folders changed.

1. **Inspect the inventory** (`landsat_scene_inventory_1980_2003.csv`). Count the months per year
   with a Landsat 5 (or L4) scene whose `clear_frac` ≥ 0.25. If years before ~1990 have long winter
   gaps, plan to fill them with climatology plus anomaly (step 4) rather than linear interpolation.
2. **Split the yearly exports**:
   `Rscript backcast/gee_1980/split_yearly_gee_exports.R <drive download folder> data/environmental_data/ndvi500_large/ndvi500_large_landsat-5`
   (and the same for `landsat-7`, `landsat-4`, `mss`, and `water_large/jrc`). File names then match
   the existing 2004+ files, so old and new files sit in the same folder.
3. **Build the per‑sensor monthly products**:
   `build_monthly_ndvi_500m("data/environmental_data/ndvi500_large/ndvi500_large_landsat-5", template, "data/environmental_data/NDVI_large/landsat-5")`
   (`environmental_data/functions/build_monthly_ndvi_500m.R`), and the same for landsat‑7 (and
   landsat‑4 / mss). This rewrites raw_mean/raw_sd/valid_fraction for 1984–2011 in one go.
4. **Harmonise** with `backcast/harmonise_landsat_to_s2.R` (`harmonise_area <- "large"`).
   Adaptations needed:
   - The L5 → adjusted L7 fit now has 1999–2011 overlap. Prefer fitting on **1999‑06 → 2003‑05
     (SLC‑on)**, or at least compare that against the 2004–2011 fit: SLC‑off stripes change which
     sub‑pixels make up the SD. Keep a hold‑out (e.g. 2002–2003, or 2010–2011 as now).
   - Add `landsat-4` (fit on adjusted L5, overlap 1984–1993) and `mss` (fit on adjusted L5/L4
     using same‑day MSS/TM pairs, 1984–1992) to the chain. Use the cover order L8 → L5 → L4 → L7 →
     MSS before 2015‑07.
   - **Gap filling.** `fill_ndvi_gaps()` (`environmental_data/functions/build_monthly_ndvi.R`)
     interpolates linearly per cell, and with `rule = 2` it holds the first value constant before
     the series starts. That is fine for 1–2 month gaps but wrong for multi‑month winter gaps in the
     1980s and for anything before 1984‑03. For gaps longer than ~2 months, use the per‑cell
     calendar‑month climatology of the harmonised series (e.g. 1985–1995) plus an anomaly. The
     anomaly comes from the nearest observed months, or from the AVHRR CDR anomaly of the
     surrounding 5 km pixel (export 04), regressed on the harmonised Landsat anomaly in 1984–2011.
     Keep `filled.tif` as the flag.
   - Fix the series start: `backcast/functions/backcast_helpers.R::fill_leading_climatology()`
     already fills months before the first layer with the mean of the first years. That is
     acceptable for the spin‑up year only.
5. **Climate**:
   - Temperature: `download_and_krig_temp(template_path = "data/environmental_data/template_raster_500_large.tif", ..., kriged_dir = "data/backcast/CDS/kriged_temperature", skip_dirs = "data/environmental_data/CDS/kriged_temperature_v2", start_date = as.Date("1979-11-01"))`
     (`environmental_data/functions/Download_and_Krig_temperature.R`). That is about 300 months of
     kriging, on 3 cores per month, so check memory (see the local machine limits). Kriging onto the
     large template also removes the 2 km western gap that `fill_edges()` now patches. If the GEE
     export (03) is used instead, write each month as the `CDSTemp_Raw` input and krig in the same
     way. Do not use raw 11 km ERA5 as temperature directly: calibration used the kriged field.
   - Precipitation: `backcast/functions/download_precipitation_backcast.py --start-year 1979 --end-year 2003`
     (CDS, same units as now). If you use the GEE file instead, divide by `days_in_month` to get
     m/day, and clamp negatives to 0.
   - `rabbit_hindcast_large.R` reads only `CDS/kriged_temperature_v2` and the single
     `precipitation.nc`. Change it to read both folders, as `Workflow_backcast_data.R` does:
     `load_monthly_stack(c(file.path(ed, "CDS", "kriged_temperature_v2"), file.path("data", "backcast", "CDS", "kriged_temperature")), template, pattern = "_Kriged\\.nc$")`,
     and the same for precipitation with `data/backcast/CDS/precipitation`.
6. **Flooding**:
   - Run `backcast/functions/download_flooding_wcs_backcast.py data/backcast/LAST_maps` (from
     1984‑09). Then run `make_monthly_floodmaps(..., min_year = 1984)`
     (`backcast/functions/make_monthly_floodmaps.R`) into the folder that
     `create_large_flood_layer.R` reads (`last_dir`), or add a second `last_dir`.
   - `create_large_flood_layer.R`: set `months <- seq(as.Date("1984-01-01"), ...)` (or Nov Y−1).
     JRC stays the off‑park source up to 2016 (`jrc_until`), using the new
     `water500_jrc_1984-03 … 2003-12` files.
   - **Before 1984‑03 (JRC) / 1984‑09 (LAST)** there is no water observation, so the script would
     leave only the no‑rabbit cells flooded. For a 1980 start, use instead a per‑cell flood
     probability per calendar month from 1984‑09 → 1994. It can be split by wet/dry hydrological
     year using ERA5‑Land Sep–Mar precipitation (`environmental_data/functions/dry_wet_years.R`
     gives the z‑scores). Set flood = 1 where the probability ≥ 0.5.
7. **No‑rabbit mask (optional)**: convert `corine1990_100m_large.tif` to the 1–44 index, then rerun
   the `rabbit_no_habitat` part of `create_large_landcover_layers.R` for 1990, or for 2000 for the
   1990s. Then either switch masks by period in `create_large_flood_layer.R` (the `blocked` cells),
   or keep 2018 and state this as a caveat.
8. **Run**: `Rscript backcast/rabbit_hindcast_large.R <n_draws> 1985 2025`, with
   `rhd_settings_list()` set to `RHD1_arrival_year = 1989` (month and seed cells to confirm with the
   park) and `RHD2_arrival_year = 2011, RHD2_arrival_month = 10`. The current script assumes RHD is
   present from the start.

## (e) Gaps, risks and recommendation

**Data gaps and risks**

- **Landsat in the 1980s and early 1990s.** Coverage of Europe/Mediterranean in the USGS archive is
  sparse and uneven before the mid‑1990s. ESA‑station scenes were repatriated, but gaps remain. Some
  months will have no scene, and winters are cloudy. The inventory (00) quantifies this before
  anything heavy is run. Landsat 5 Tier 2 and Landsat 4 can add scenes. Linear gap filling across
  long gaps biases the winter/spring NDVI that drives K in the breeding season.
- **Harmonisation chain.** Landsat 5 reaches Sentinel‑2 only through Landsat 7; L4 and MSS add
  another step each. Errors add up and are not in the posterior spread. Report the hold‑out RMSE per
  step (`harmonisation_fits.csv`), as for phase 2. MSS SD (60 m pixels, TOA) is not the S2 10 m SD:
  treat SD before 1984 as essentially climatology.
- **AVHRR** sensor changes (NOAA‑7/9/11/14, orbital drift) cause artefacts at 5 km. Use it only as
  a gap‑filling anomaly, never as the level.
- **Flooding.** LAST starts on 1984‑09‑10 and JRC on 1984‑03. 1980–83 flooding can only be modelled.
  JRC misses about 40 % of flooded cells relative to LAST (kappa 0.63, from the 2004–2025
  comparison), probably worse with Landsat 5 alone.
- **Static layers.** Vegetation height is a ~2015–2020 map. The 1980s–90s saw large shrub and
  forest changes: eucalyptus removal in El Abalario, pine plantations, and fires such as Moguer 2017
  (which falls after 2005). Outside the park, greenhouses and berry farms and built‑up areas
  expanded strongly after 1990, so a CLC 2018 mask understates habitat in the 1980s. CLC 1990/2000
  partly fixes this. There is no historical height product.
- **Climate.** ERA5‑Land is one consistent product, so no harmonisation is needed. But 1980s
  station density in the reanalysis is lower, and the kriging only adds DEM‑driven detail, not new
  information.
- **Disease and extrapolation.** The posterior is from 2017–2024 (post‑RHD2) park census data. Before
  RHD1 (~1988–89) rabbits in Doñana were several times more abundant: they are now under 10 % of
  pre‑RHD numbers. The model must reproduce that with the K map plus the RHD module, and
  myxomatosis (endemic since the 1950s) is not modelled. Results before 2005 (census start) cannot
  be validated. Pre‑RHD densities depend on the RHD1 arrival date and mortality settings, which
  still need calibration (`backcast_report.md`, RHD section).

**Recommendation: 1985 (covering 1990), not 1980.**

| | Start 1990 | Start 1985 (recommended) | Start 1980 |
|---|---|---|---|
| NDVI | Landsat 5 TM, decent coverage | Landsat 5 TM from 1984‑03 (spin‑up 1985 absorbs sparse 1984) | MSS/AVHRR only for 1980–84‑02; nothing at all before 1981‑06 except MSS |
| Flooding | LAST + JRC | LAST + JRC (from 1984‑09 / 03) | modelled climatology 1980–84 |
| Climate | ERA5‑Land | ERA5‑Land | ERA5‑Land |
| RHD | arrival just before the start (needs an initial immune share) | pre‑RHD baseline 1986–88, then the RHD1 wave | longer pre‑RHD period, but on weak inputs |
| GEE tasks | ~48 (same set; Landsat could start 1988) | ~48 | ~77 |

Starting in 1985 costs nothing extra compared with 1990, because the same exports serve both. It
also includes the pre‑RHD years, which give the RHD1 decline target (about 60 %). A 1980 start adds
four years whose NDVI and flooding are largely reconstructed rather than observed. It is worth doing
only as a sensitivity run once the 1985 run works, and only if the inventory shows usable MSS
coverage for 1980–83.

**Sources** (checked 8 Nov 2026): GEE catalogue pages for
[ECMWF/ERA5_LAND/MONTHLY_AGGR](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY_AGGR),
[LANDSAT/LT05/C02/T1_L2](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C02_T1_L2),
[LANDSAT/LT04/C02/T1_L2](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C02_T1_L2),
[LM02](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LM02_C02_T1),
[LM03](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LM03_C02_T1),
[LM04](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LM04_C02_T1),
[LM05](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LM05_C02_T1),
[NOAA/CDR/AVHRR/NDVI/V5](https://developers.google.com/earth-engine/datasets/catalog/NOAA_CDR_AVHRR_NDVI_V5),
[NASA/GIMMS/3GV0](https://developers.google.com/earth-engine/datasets/catalog/NASA_GIMMS_3GV0),
[JRC/GSW1_4/MonthlyHistory](https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_MonthlyHistory),
[IDAHO_EPSCOR/TERRACLIMATE](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE),
[COPERNICUS/CORINE/V20/100m](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_CORINE_V20_100m).
On Landsat coverage of Europe in the 1980s and 1990s:
[Goward et al., historical record of Landsat global coverage](https://webarchive.library.unt.edu/eot2008/20080916135403/http://neptune.gsfc.nasa.gov/publications/pdf/pubs2006/Goward-Historical_Record_of_Landsat_Global_Coverage.pdf),
[Landsat/Sentinel‑2 availability over Europe (preprint)](https://www.preprints.org/frontend/manuscript/c9fa3c63932e329c960b4dc3ad8ad493/download_pub).
