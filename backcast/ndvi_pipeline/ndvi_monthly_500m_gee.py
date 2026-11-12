"""
Monthly 500 m NDVI statistics computed on Google Earth Engine (same products as ndvi_monthly_500m.py).

Per month and sensor, on GEE:
  1. all scenes of the month (scene cloud cover < --cloud-cover), cloud-masked per pixel
     (Sentinel-2: SCL classes 2, 4, 5, 6, 7; Landsat: QA_PIXEL bits 0, 1, 3, 4, 5), NDVI per scene
  2. monthly composite: maximum NDVI per 10 m pixel (Landsat 30 m pixels kept intact, nearest)
  3. per 500 m template cell: mean, SD (n - 1) and valid fraction (share of the 50 x 50 pixels
     with a value)
and exports one GeoTIFF per month to a Google Drive folder:
  ndvi500_<sensor>_<YYYY-MM>.tif     bands: mean, sd, valid_fraction; -9999 = no data
Months without any scene are not exported (they become gaps, filled in R).

The 10 m grid is the template grid with its origin rounded to whole 10 m, as in
ndvi_monthly_500m.py, in WGS 84 / UTM 29N. The exported grid is therefore the template grid shifted
by < 10 m; build_monthly_ndvi_500m.R puts it back on the template.

Sentinel-2 uses COPERNICUS/S2_SR_HARMONIZED, in which the reflectance offset of processing
baseline 04.00 (from Jan 2022) is already removed. Landsat uses Collection 2 Tier 1 Level-2.

Setup (once): Earth Engine account and Cloud project, `pip install earthengine-api`,
`earthengine authenticate`. Then, from the Donana_Digital_Twin project root, e.g.:
  python backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py --project my-ee-project \
      --template data/environmental_data/template_raster_500_large.tif \
      --sensor sentinel-2 --start 2015-07 --end 2026-10 --drive-folder ndvi500_large_sentinel-2
This only submits export tasks (the computation runs on GEE). Follow them at
https://code.earthengine.google.com/tasks or with `earthengine task list`. Months whose file is
already in --skip-dir (a local folder with downloaded files) are not submitted again.
"""
import argparse
import calendar
import os
import sys

import ee
import rasterio

FINE_RES = 10
NODATA = -9999
OUTPUT_CRS = "EPSG:32629"   # WGS 84 / UTM 29N, as the pipeline; the template is ETRS89 / UTM 29N (< 1 m apart)

S2_VALID_SCL = [2, 4, 5, 6, 7]
LANDSAT_BAD_BITS = [0, 1, 3, 4, 5]   # fill, dilated cloud, cloud, cloud shadow, snow
LANDSAT_SCALE, LANDSAT_OFFSET = 0.0000275, -0.2

SENSORS = {
    "sentinel-2": dict(collection="COPERNICUS/S2_SR_HARMONIZED", cloud="CLOUDY_PIXEL_PERCENTAGE"),
    "landsat-5": dict(collection="LANDSAT/LT05/C02/T1_L2", cloud="CLOUD_COVER", red="SR_B3", nir="SR_B4"),
    "landsat-7": dict(collection="LANDSAT/LE07/C02/T1_L2", cloud="CLOUD_COVER", red="SR_B3", nir="SR_B4"),
    "landsat-8": dict(collection="LANDSAT/LC08/C02/T1_L2", cloud="CLOUD_COVER", red="SR_B4", nir="SR_B5"),
    "landsat-9": dict(collection="LANDSAT/LC09/C02/T1_L2", cloud="CLOUD_COVER", red="SR_B4", nir="SR_B5"),
}


def template_grid(path):
    with rasterio.open(path) as src:
        t = src.transform
        if abs(t.e + t.a) > 1e-6 or t.b or t.d:
            raise ValueError("Template must be north-up with square cells")
        fact = int(round(t.a / FINE_RES))
        if fact * FINE_RES != t.a:
            raise ValueError(f"Template resolution {t.a} is not a multiple of {FINE_RES} m")
        x0 = FINE_RES * round(t.c / FINE_RES)
        y0 = FINE_RES * round(t.f / FINE_RES)
        return dict(x0=x0, y0=y0, res=t.a, fact=fact, ncol=src.width, nrow=src.height)


def s2_ndvi(img):
    scl = img.select("SCL")
    valid = scl.remap(S2_VALID_SCL, [1] * len(S2_VALID_SCL), 0).eq(1)
    red = img.select("B4").toFloat()
    nir = img.select("B8").toFloat()
    den = nir.add(red)
    ndvi = nir.subtract(red).divide(den).updateMask(den.neq(0)).clamp(-1, 1)
    return ndvi.updateMask(valid).rename("ndvi")


def landsat_ndvi(spec):
    bad = sum(1 << b for b in LANDSAT_BAD_BITS)

    def f(img):
        valid = img.select("QA_PIXEL").bitwiseAnd(bad).eq(0)
        red = img.select(spec["red"]).multiply(LANDSAT_SCALE).add(LANDSAT_OFFSET)
        nir = img.select(spec["nir"]).multiply(LANDSAT_SCALE).add(LANDSAT_OFFSET)
        den = nir.add(red)
        ndvi = nir.subtract(red).divide(den).updateMask(den.neq(0)).clamp(-1, 1)
        return ndvi.updateMask(valid).rename("ndvi")
    return f


def month_image(sensor, year, month, region, grid, cloud_cover):
    spec = SENSORS[sensor]
    start = ee.Date.fromYMD(year, month, 1)
    col = (ee.ImageCollection(spec["collection"])
           .filterBounds(region)
           .filterDate(start, start.advance(1, "month")))
    if cloud_cover is not None:
        col = col.filter(ee.Filter.lt(spec["cloud"], cloud_cover))
    n = col.size()
    ndvi = col.map(s2_ndvi if sensor == "sentinel-2" else landsat_ndvi(spec))

    fine = [FINE_RES, 0, grid["x0"], 0, -FINE_RES, grid["y0"]]
    coarse = [grid["res"], 0, grid["x0"], 0, -grid["res"], grid["y0"]]
    # The composite gets the 10 m grid as its default projection, so reduceResolution aggregates
    # 10 m pixels (nearest from the source). Not .reproject(): a forced reprojection of the whole
    # area at 10 m is larger than GEE allows ("Reprojection output too large").
    comp = ndvi.max().setDefaultProjection(crs=OUTPUT_CRS, crsTransform=fine)

    reducer = (ee.Reducer.mean()
               .combine(ee.Reducer.stdDev(), sharedInputs=True)
               .combine(ee.Reducer.count(), sharedInputs=True))
    stats = comp.reduceResolution(reducer=reducer, maxPixels=grid["fact"] ** 2 + 100) \
                .reproject(crs=OUTPUT_CRS, crsTransform=coarse)
    count = stats.select("ndvi_count")
    # GEE's stdDev is the population SD; R's sd() (used so far) is the sample SD
    sd = stats.select("ndvi_stdDev").multiply(count.divide(count.subtract(1)).sqrt()).updateMask(count.gt(1))
    valid = count.divide(grid["fact"] ** 2).unmask(0)
    out = ee.Image.cat([stats.select("ndvi_mean").rename("mean"), sd.rename("sd"),
                        valid.rename("valid_fraction")]).toFloat().unmask(NODATA)
    return out, n


def months_between(start, end):
    y, m = map(int, start.split("-"))
    y1, m1 = map(int, end.split("-"))
    while (y, m) <= (y1, m1):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--project", required=True, help="Google Cloud project registered for Earth Engine")
    p.add_argument("--template", required=True)
    p.add_argument("--sensor", required=True, choices=sorted(SENSORS))
    p.add_argument("--start", required=True, metavar="YYYY-MM")
    p.add_argument("--end", required=True, metavar="YYYY-MM")
    p.add_argument("--drive-folder", required=True, help="Google Drive folder for the exports")
    p.add_argument("--cloud-cover", type=float, default=60, help="Scene cloud cover limit in %% (default 60; -1 = none)")
    p.add_argument("--shard-cells", type=int, default=40,
                   help="Export shard size in template cells (default 40 = 2000 x 2000 pixels at 10 m)")
    p.add_argument("--skip-dir", default=None, help="Local folder with already downloaded files: those months are skipped")
    a = p.parse_args()

    ee.Initialize(project=a.project)
    grid = template_grid(a.template)
    proj = ee.Projection(OUTPUT_CRS)
    x1 = grid["x0"] + grid["ncol"] * grid["res"]
    y1 = grid["y0"] - grid["nrow"] * grid["res"]
    region = ee.Geometry.Rectangle([grid["x0"], y1, x1, grid["y0"]], proj, False)
    cloud = None if a.cloud_cover < 0 else a.cloud_cover

    submitted, empty, skipped = 0, 0, 0
    for y, m in months_between(a.start, a.end):
        name = f"ndvi500_{a.sensor}_{y:04d}-{m:02d}"
        if a.skip_dir and os.path.exists(os.path.join(a.skip_dir, name + ".tif")):
            skipped += 1
            continue
        img, n = month_image(a.sensor, y, m, region, grid, cloud)
        n = n.getInfo()
        if n == 0:
            print(f"{name}: no scenes, not exported")
            empty += 1
            continue
        task = ee.batch.Export.image.toDrive(
            image=img, description=name, folder=a.drive_folder, fileNamePrefix=name,
            crs=OUTPUT_CRS, crsTransform=[grid["res"], 0, grid["x0"], 0, -grid["res"], grid["y0"]],
            region=region, maxPixels=1e9, fileFormat="GeoTIFF",
            # GEE computes the export in shards of shardSize x shardSize output cells, each needing
            # (shardSize * fact)^2 pixels at 10 m. The default (256) asks for the whole area at
            # 10 m at once: "Reprojection output too large".
            shardSize=a.shard_cells)
        task.start()
        submitted += 1
        print(f"{name}: {n} scenes, task submitted")
    print(f"Submitted {submitted}, no scenes {empty}, already downloaded {skipped}")


if __name__ == "__main__":
    sys.exit(main())
