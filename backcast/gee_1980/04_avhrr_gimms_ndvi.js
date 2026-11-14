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
