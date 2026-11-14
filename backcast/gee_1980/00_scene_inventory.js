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
