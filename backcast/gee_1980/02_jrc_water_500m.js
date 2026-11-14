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
