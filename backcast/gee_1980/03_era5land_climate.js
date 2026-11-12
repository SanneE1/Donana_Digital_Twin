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
