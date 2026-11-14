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
