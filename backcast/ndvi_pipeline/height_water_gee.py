"""
Vegetation height and monthly surface water for a large template, computed on Google Earth Engine.

The rabbit model needs vegetation height (static) and flooding (monthly) in every cell. The park
maps (LiDAR vegetation height, LAST flood maps) only cover the park, so for the large (lynx) template
these come from global products, checked against the park maps in R.

Two commands, both exporting GeoTIFFs to a Google Drive folder on the template grid (origin rounded
to whole 10 m, WGS 84 / UTM 29N, as ndvi_monthly_500m_gee.py; -9999 = no data):

  height   one file per product, on a 100 m grid (5 x 5 per template cell; aggregated in R):
             height100_<product>.tif   band: mean canopy height in m
           products: eth  (Lang et al. 2023, 10 m, 2020)
                     meta (Tolan et al. 2024, 1 m, ~2018-2020)
                     glad (Potapov et al. 2021, 30 m, 2019)

  water    one file per month and source, on the template grid:
             water500_<source>_<YYYY-MM>.tif   bands: water_fraction (share of observed pixels that
                                               were water at least once in the month),
                                               valid_fraction (share of pixels observed)
           sources: jrc          JRC Global Surface Water monthly history (30 m, to 2021)
                    sentinel-2   SCL class 6 (water) among clear pixels (SCL 2, 4, 5, 6, 7)
                    landsat-5/7/8  QA_PIXEL bit 7 (water) among clear pixels

From the Donana_Digital_Twin project root, e.g.:
  python backcast/ndvi_pipeline/height_water_gee.py height --project my-ee-project \
      --template data/environmental_data/template_raster_500_large.tif --drive-folder height_large
  python backcast/ndvi_pipeline/height_water_gee.py water --project my-ee-project \
      --template data/environmental_data/template_raster_500_large.tif \
      --source jrc --start 2004-01 --end 2021-12 --drive-folder water500_large_jrc
"""
import argparse
import os
import sys

import ee

from ndvi_monthly_500m_gee import (FINE_RES, NODATA, OUTPUT_CRS, S2_VALID_SCL, LANDSAT_BAD_BITS,
                                   SENSORS, months_between, template_grid)

def _mosaic(collection_id):
    # a mosaic gets a 1-degree default projection; reduceResolution needs the native pixels
    col = ee.ImageCollection(collection_id)
    return col.mosaic().setDefaultProjection(col.first().projection())


HEIGHT_PRODUCTS = {
    "eth":  lambda: ee.Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1").select(0),
    "meta": lambda: _mosaic("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").select(0),
    "glad": lambda: ee.Image("users/potapovpeter/GEDI_V27").select(0),
}
# native resolution in m, for the reduceResolution pixel limit
HEIGHT_NATIVE_RES = {"eth": 10, "meta": 1, "glad": 30}
HEIGHT_RES = 100
LANDSAT_WATER_BIT = 7


def region_of(grid):
    proj = ee.Projection(OUTPUT_CRS)
    x1 = grid["x0"] + grid["ncol"] * grid["res"]
    y1 = grid["y0"] - grid["nrow"] * grid["res"]
    return ee.Geometry.Rectangle([grid["x0"], y1, x1, grid["y0"]], proj, False)


def export(img, name, folder, grid, res, region, shard):
    task = ee.batch.Export.image.toDrive(
        image=img, description=name, folder=folder, fileNamePrefix=name,
        crs=OUTPUT_CRS, crsTransform=[res, 0, grid["x0"], 0, -res, grid["y0"]],
        region=region, maxPixels=1e9, fileFormat="GeoTIFF", shardSize=shard)
    task.start()


def cmd_height(a, grid, region):
    for prod in a.products:
        img = HEIGHT_PRODUCTS[prod]()
        n_px = int((HEIGHT_RES / HEIGHT_NATIVE_RES[prod]) ** 2) + 100
        mean = (img.toFloat()
                .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=n_px)
                .reproject(crs=OUTPUT_CRS, crsTransform=[HEIGHT_RES, 0, grid["x0"], 0, -HEIGHT_RES, grid["y0"]])
                .rename("height").unmask(NODATA).toFloat())
        name = f"height100_{prod}"
        export(mean, name, a.drive_folder, grid, HEIGHT_RES, region, 64)
        print(f"{name}: task submitted")


def water_month(source, year, month, region, grid, cloud_cover):
    """Per fine pixel: 1 if water in any clear observation of the month, 0 if observed and never water."""
    start = ee.Date.fromYMD(year, month, 1)
    end = start.advance(1, "month")
    if source == "jrc":
        col = ee.ImageCollection("JRC/GSW1_4/MonthlyHistory").filterDate(start, end)
        n = col.size()
        w = col.first().select("water")
        water = w.eq(2).updateMask(w.gt(0))
        native = 30
    else:
        spec = SENSORS[source]
        col = ee.ImageCollection(spec["collection"]).filterBounds(region).filterDate(start, end)
        if cloud_cover is not None:
            col = col.filter(ee.Filter.lt(spec["cloud"], cloud_cover))
        n = col.size()
        if source == "sentinel-2":
            def f(img):
                scl = img.select("SCL")
                valid = scl.remap(S2_VALID_SCL, [1] * len(S2_VALID_SCL), 0).eq(1)
                return scl.eq(6).updateMask(valid).rename("water")
        else:
            bad = sum(1 << b for b in LANDSAT_BAD_BITS)

            def f(img):
                qa = img.select("QA_PIXEL")
                valid = qa.bitwiseAnd(bad).eq(0)
                return qa.bitwiseAnd(1 << LANDSAT_WATER_BIT).neq(0).updateMask(valid).rename("water")
        water = col.map(f).max()
        native = 10
    fine = [FINE_RES, 0, grid["x0"], 0, -FINE_RES, grid["y0"]]
    coarse = [grid["res"], 0, grid["x0"], 0, -grid["res"], grid["y0"]]
    # as in ndvi_monthly_500m_gee.py: 10 m default projection, nearest from the source pixels
    water = water.rename("water").toFloat().setDefaultProjection(crs=OUTPUT_CRS, crsTransform=fine)
    stats = water.reduceResolution(
        reducer=ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True),
        maxPixels=grid["fact"] ** 2 + 100).reproject(crs=OUTPUT_CRS, crsTransform=coarse)
    valid = stats.select("water_count").divide(grid["fact"] ** 2).unmask(0)
    out = ee.Image.cat([stats.select("water_mean").rename("water_fraction"),
                        valid.rename("valid_fraction")]).toFloat().unmask(NODATA)
    return out, n


def cmd_water(a, grid, region):
    cloud = None if a.cloud_cover < 0 else a.cloud_cover
    submitted, empty, skipped = 0, 0, 0
    for y, m in months_between(a.start, a.end):
        name = f"water500_{a.source}_{y:04d}-{m:02d}"
        if a.skip_dir and os.path.exists(os.path.join(a.skip_dir, name + ".tif")):
            skipped += 1
            continue
        img, n = water_month(a.source, y, m, region, grid, cloud)
        n = n.getInfo()
        if n == 0:
            print(f"{name}: no scenes, not exported")
            empty += 1
            continue
        export(img, name, a.drive_folder, grid, grid["res"], region, a.shard_cells)
        submitted += 1
        print(f"{name}: {n} scenes, task submitted")
    print(f"Submitted {submitted}, no scenes {empty}, already downloaded {skipped}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("height", "water"):
        s = sub.add_parser(name)
        s.add_argument("--project", required=True, help="Google Cloud project registered for Earth Engine")
        s.add_argument("--template", required=True)
        s.add_argument("--drive-folder", required=True, help="Google Drive folder for the exports")
    sub.choices["height"].add_argument("--products", nargs="+", default=sorted(HEIGHT_PRODUCTS),
                                       choices=sorted(HEIGHT_PRODUCTS))
    w = sub.choices["water"]
    w.add_argument("--source", required=True, choices=["jrc", "sentinel-2", "landsat-5", "landsat-7", "landsat-8"])
    w.add_argument("--start", required=True, metavar="YYYY-MM")
    w.add_argument("--end", required=True, metavar="YYYY-MM")
    w.add_argument("--cloud-cover", type=float, default=60, help="Scene cloud cover limit in %% (default 60; -1 = none)")
    w.add_argument("--shard-cells", type=int, default=40)
    w.add_argument("--skip-dir", default=None, help="Local folder with already downloaded files: those months are skipped")
    a = p.parse_args()

    ee.Initialize(project=a.project)
    grid = template_grid(a.template)
    region = region_of(grid)
    return cmd_height(a, grid, region) if a.cmd == "height" else cmd_water(a, grid, region)


if __name__ == "__main__":
    sys.exit(main())
