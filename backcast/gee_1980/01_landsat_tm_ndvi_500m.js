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
