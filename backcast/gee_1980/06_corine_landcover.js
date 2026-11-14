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
