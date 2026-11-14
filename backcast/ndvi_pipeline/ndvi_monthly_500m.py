"""
Monthly 500 m NDVI statistics straight from Sentinel-2 or Landsat, one month and one tile at a time.

Produces, per month, the same raw products as build_monthly_ndvi.R makes from the 10 m files of
pipeline.py, without writing the 10 m data:
  1. all scenes of the month on a 10 m grid aligned with the template (as the "near" resampling
     in build_monthly_ndvi.R), cloud-masked, NDVI per scene
  2. monthly composite: the maximum NDVI per 10 m pixel
  3. per 500 m template cell: mean and SD (n - 1) of the composite pixels, and the valid fraction
     (share of the 50 x 50 pixels with a value)

Memory use depends on the tile size (--tile-cells), not on the size of the area or the number of
scenes, and the only files written are small GeoTIFFs on the template grid:
  <output-dir>/ndvi500_<sensor>_<YYYY-MM>.tif    bands: mean, sd, valid_fraction
  <output-dir>/ndvi500_<sensor>_<YYYY-MM>.none   month without any scene (not tried again)
Months already on disk are skipped, so the script can be stopped and restarted. A month in which a
tile fails after several attempts is not written, and is tried again on the next run.

Gap filling and stacking into the model inputs: environmental_data/functions/build_monthly_ndvi_500m.R

Examples (from the Donana_Digital_Twin project root):
  python backcast/ndvi_pipeline/ndvi_monthly_500m.py --template data/environmental_data/template_raster_500_large.tif \
      --source mpc --start 2015-07 --end 2026-10 --output-dir data/environmental_data/ndvi500_large/sentinel-2
  python backcast/ndvi_pipeline/ndvi_monthly_500m.py --template ... --source landsat --platforms landsat-7 \
      --start 2005-01 --end 2015-06 --output-dir data/environmental_data/ndvi500_large/landsat-7
"""
import argparse
import calendar
import json
import logging
import math
import os
import sys
import time
import warnings

import _proj_env  # noqa: F401 — must precede other geospatial imports

import dask
import numpy as np
import planetary_computer
import pystac_client
import rasterio
from affine import Affine
from odc.geo.geobox import GeoBox
from odc.stac import load as odc_load
from pyproj import Transformer

import config
import pipeline as pl

log = logging.getLogger("ndvi500")

S2_BOA_ADD_OFFSET = 1000   # processing baseline >= 04.00: DN = reflectance * 10000 + 1000
FINE_RES = 10              # m, grid of the monthly composite (as RESOLUTION in config.py)


def has_boa_offset(item) -> bool:
    """True when a Sentinel-2 item's reflectances still carry the baseline 04.00 offset."""
    if item.properties.get("earthsearch:boa_offset_applied"):
        return False
    baseline = item.properties.get("s2:processing_baseline") or item.properties.get("processing:version")
    if baseline is None:
        raise ValueError(f"No processing baseline for item {item.id}: cannot tell whether to remove the offset")
    return float(baseline) >= 4.0


# ── Grid ─────────────────────────────────────────────────────────────────────

class Template:
    def __init__(self, path):
        with rasterio.open(path) as src:
            self.transform = src.transform
            self.crs = src.crs
            self.width, self.height = src.width, src.height
            self.land = ~np.isnan(src.read(1, masked=True).filled(np.nan).astype("float64"))
        self.res = self.transform.a
        if abs(self.transform.e + self.res) > 1e-6 or self.transform.b or self.transform.d:
            raise ValueError("Template must be north-up with square cells")
        self.fact = int(round(self.res / FINE_RES))
        if self.fact * FINE_RES != self.res:
            raise ValueError(f"Template resolution {self.res} is not a multiple of {FINE_RES} m")

    def tiles(self, size):
        for r0 in range(0, self.height, size):
            for c0 in range(0, self.width, size):
                nr, nc = min(size, self.height - r0), min(size, self.width - c0)
                yield r0, c0, nr, nc, self.land[r0:r0 + nr, c0:c0 + nc].any()

    def fine_geobox(self, r0, c0, nr, nc):
        """10 m grid of a tile. The origin is rounded to whole 10 m, so the source pixels map 1:1 to
        the template-aligned 10 m pixels, as with the nearest-neighbour resampling in
        build_monthly_ndvi.R (the template is offset from the 10 m source grid by < 5 m)."""
        x0 = self.transform.c + c0 * self.res
        y0 = self.transform.f - r0 * self.res
        x0 = FINE_RES * round(x0 / FINE_RES)
        y0 = FINE_RES * round(y0 / FINE_RES)
        aff = Affine(FINE_RES, 0, x0, 0, -FINE_RES, y0)
        # Data are loaded in the pipeline's CRS (WGS 84 / UTM 29N); the template is ETRS89 / UTM 29N.
        # They differ by < 1 m here and build_monthly_ndvi.R treats them as the same, as done here.
        return GeoBox((nr * self.fact, nc * self.fact), aff, config.OUTPUT_CRS)

    def bbox_wgs84(self, r0=0, c0=0, nr=None, nc=None, margin=0.0):
        nr = self.height if nr is None else nr
        nc = self.width if nc is None else nc
        x0 = self.transform.c + c0 * self.res - margin
        x1 = self.transform.c + (c0 + nc) * self.res + margin
        y1 = self.transform.f - r0 * self.res + margin
        y0 = self.transform.f - (r0 + nr) * self.res - margin
        tr = Transformer.from_crs(self.crs, "EPSG:4326", always_xy=True)
        xs, ys = tr.transform([x0, x0, x1, x1], [y0, y1, y0, y1])
        return [min(xs), min(ys), max(xs), max(ys)]


def overlaps(item, bbox):
    w, s, e, n = item.bbox
    return not (e < bbox[0] or w > bbox[2] or n < bbox[1] or s > bbox[3])


# ── One tile ─────────────────────────────────────────────────────────────────

def _load(items, geobox, bands, resampling):
    return odc_load(items, geobox=geobox, bands=bands, groupby="solar_day",
                    chunks={"time": 1, "x": 1024, "y": 1024}, dtype="float32", nodata=0,
                    resampling=resampling)


def _ndvi(red, nir):
    den = nir + red
    return ((nir - red) / den.where(den != 0)).clip(-1, 1)


def composite_sentinel2(items, geobox):
    """Maximum NDVI per pixel over all scenes. Scenes with and without the reflectance offset are
    loaded separately (one day can mix processing baselines), as in NDVI-Pipeline-Donana-main_v2."""
    groups = {}
    for it in items:
        groups.setdefault(has_boa_offset(it), []).append(it)
    out = None
    for offset, group in groups.items():
        ds = _load(group, geobox, config.BANDS, {"SCL": "nearest", "*": "average"})
        ds = ds.where(ds != 0)
        red, nir = ds["B04"], ds["B08"]
        if offset:
            red = (red - S2_BOA_ADD_OFFSET).clip(min=0)
            nir = (nir - S2_BOA_ADD_OFFSET).clip(min=0)
        valid = ds["SCL"].isin(config.VALID_SCL_CLASSES)
        comp = _ndvi(red.where(valid), nir.where(valid)).max("time", skipna=True).values
        out = comp if out is None else np.fmax(out, comp)
    return out


def composite_landsat(items, geobox):
    """Maximum NDVI per pixel over all scenes; 30 m pixels kept intact on the 10 m grid (nearest)."""
    ds = _load(items, geobox, config.LANDSAT_BANDS, "nearest")
    qa = ds["QA"].fillna(1).astype("uint16")
    bad = sum(1 << b for b in config.LANDSAT_QA_MASK_BITS)
    valid = (qa & bad) == 0
    red = (ds["B04"] * config.LANDSAT_SCALE + config.LANDSAT_OFFSET).where(valid)
    nir = (ds["B08"] * config.LANDSAT_SCALE + config.LANDSAT_OFFSET).where(valid)
    return _ndvi(red, nir).max("time", skipna=True).values


def block_stats(comp, fact):
    """Mean, SD (n - 1) and valid fraction of each fact x fact block."""
    nr, nc = comp.shape[0] // fact, comp.shape[1] // fact
    blocks = comp.reshape(nr, fact, nc, fact).transpose(0, 2, 1, 3).reshape(nr, nc, fact * fact)
    count = np.sum(~np.isnan(blocks), axis=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(blocks, axis=2)
        sd = np.nanstd(blocks, axis=2, ddof=1)
    sd[count < 2] = np.nan
    return mean, sd, count / (fact * fact)


# ── One month ────────────────────────────────────────────────────────────────

def month_window(year, month):
    last = calendar.monthrange(year, month)[1]
    return f"{year:04d}-{month:02d}-01T00:00:00Z/{year:04d}-{month:02d}-{last:02d}T23:59:59Z"


def sensor_name():
    if config.DATA_SOURCE == "landsat":
        return "+".join(config.LANDSAT_PLATFORMS)
    return "sentinel-2"


def process_month(tpl, year, month, out_dir, tile_cells, attempts):
    tag = f"{year:04d}-{month:02d}"
    base = os.path.join(out_dir, f"ndvi500_{sensor_name()}_{tag}")
    if os.path.exists(base + ".tif") or os.path.exists(base + ".none"):
        return "skipped"

    config.DATE_RANGE = month_window(year, month)
    # Unsigned search: the items are signed per tile attempt, as the signed URLs expire
    client = pystac_client.Client.open(config.MPC_STAC_URL)
    bbox = tpl.bbox_wgs84()
    query = {}
    if config.MAX_CLOUD_COVER is not None:
        query["eo:cloud_cover"] = {"lt": config.MAX_CLOUD_COVER}
    if config.DATA_SOURCE == "landsat":
        query["platform"] = {"in": config.LANDSAT_PLATFORMS}
    items = list(client.search(collections=[pl.get_collection_name()], bbox=bbox,
                               datetime=config.DATE_RANGE, query=query or None).items())
    if not items:
        with open(base + ".none", "w") as f:
            json.dump({"month": tag, "sensor": sensor_name(), "reason": "no scenes"}, f)
        log.info(f"{tag}: no scenes")
        return "none"
    items = pl.normalise_band_names(pl._ensure_proj_extension(items))
    days = sorted({it.datetime.date() for it in items})
    log.info(f"{tag}: {len(items)} items on {len(days)} days")

    shape = (tpl.height, tpl.width)
    mean, sd, valid = (np.full(shape, np.nan, dtype="float32") for _ in range(3))
    composite = composite_landsat if config.DATA_SOURCE == "landsat" else composite_sentinel2

    t0 = time.time()
    for r0, c0, nr, nc, has_land in tpl.tiles(tile_cells):
        if not has_land:
            continue
        tile_items = [it for it in items if overlaps(it, tpl.bbox_wgs84(r0, c0, nr, nc, margin=100))]
        if not tile_items:
            continue
        gb = tpl.fine_geobox(r0, c0, nr, nc)
        for attempt in range(1, attempts + 1):
            try:
                comp = composite([planetary_computer.sign(it) for it in tile_items], gb)
                break
            except Exception as e:  # network errors, expired tokens, throttling
                log.warning(f"{tag} tile ({r0},{c0}) attempt {attempt} failed: {e}")
                if attempt == attempts:
                    log.error(f"{tag}: giving up on this month, it is tried again on the next run")
                    return "failed"
                time.sleep(20 * attempt)
        m, s, v = block_stats(comp, tpl.fact)
        mean[r0:r0 + nr, c0:c0 + nc] = m
        sd[r0:r0 + nr, c0:c0 + nc] = s
        valid[r0:r0 + nr, c0:c0 + nc] = v

    profile = dict(driver="GTiff", width=tpl.width, height=tpl.height, count=3, dtype="float32",
                   crs=tpl.crs, transform=tpl.transform, nodata=np.nan, compress="deflate")
    tmp = base + ".tif.part"
    with rasterio.open(tmp, "w", **profile) as dst:
        for i, (arr, name) in enumerate(zip((mean, sd, valid), ("mean", "sd", "valid_fraction")), start=1):
            dst.write(arr, i)
            dst.set_band_description(i, name)
        dst.update_tags(month=tag, sensor=sensor_name(), n_items=len(items), n_days=len(days),
                        min_valid_frac_applied="no", composite="monthly maximum NDVI per 10 m pixel",
                        reflectance_offset_removed="yes")
    os.replace(tmp, base + ".tif")
    log.info(f"{tag}: written in {time.time() - t0:.0f} s")
    return "done"


# ── CLI ──────────────────────────────────────────────────────────────────────

def months_between(start, end):
    y, m = map(int, start.split("-"))
    y1, m1 = map(int, end.split("-"))
    while (y, m) <= (y1, m1):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--template", required=True, help="500 m template raster (grid of the output)")
    p.add_argument("--start", required=True, metavar="YYYY-MM")
    p.add_argument("--end", required=True, metavar="YYYY-MM")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--source", choices=["mpc", "landsat"], default="mpc",
                   help="mpc: Sentinel-2 L2A on Planetary Computer; landsat: Landsat C2 L2 on Planetary Computer")
    p.add_argument("--platforms", nargs="+", default=None, help="Landsat platforms, e.g. landsat-7")
    p.add_argument("--cloud-cover", type=float, default=None,
                   help="Scene-level cloud cover limit in %% (default: config.py; -1 = none)")
    p.add_argument("--tile-cells", type=int, default=40,
                   help="Tile size in template cells (default 40 = 20 km, about 0.5-1 GB of memory)")
    p.add_argument("--workers", type=int, default=4, help="Dask threads (default 4)")
    p.add_argument("--attempts", type=int, default=3, help="Attempts per tile before a month is left for the next run")
    a = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", force=True)
    config.DATA_SOURCE = a.source
    if a.platforms:
        config.LANDSAT_PLATFORMS = a.platforms
    if a.cloud_cover is not None:
        config.MAX_CLOUD_COVER = None if a.cloud_cover < 0 else a.cloud_cover
    dask.config.set(scheduler="threads", num_workers=a.workers)
    os.makedirs(a.output_dir, exist_ok=True)

    tpl = Template(a.template)
    log.info(f"Template {tpl.width} x {tpl.height} cells of {tpl.res:g} m; sensor {sensor_name()}; "
             f"tiles of {a.tile_cells} cells")
    status = {}
    for y, m in months_between(a.start, a.end):
        status[f"{y:04d}-{m:02d}"] = process_month(tpl, y, m, a.output_dir, a.tile_cells, a.attempts)
    try:
        import psutil
        mi = psutil.Process().memory_info()
        peak = getattr(mi, "peak_wset", None) or getattr(mi, "rss", 0)
        log.info(f"Peak memory: {peak / 2**30:.2f} GB")
    except ImportError:
        pass
    failed = [k for k, v in status.items() if v == "failed"]
    log.info("Finished: " + ", ".join(f"{v} {sum(1 for s in status.values() if s == v)}"
                                       for v in ("done", "skipped", "none", "failed")))
    if failed:
        log.error("Failed months (run again to retry): " + ", ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
