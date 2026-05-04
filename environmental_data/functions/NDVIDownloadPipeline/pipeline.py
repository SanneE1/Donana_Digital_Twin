"""
Sentinel-2 NDVI data cube pipeline.

AOI extraction → STAC search → load → cloud mask → NDVI → export.
"""
import argparse
import os
import pyproj

print(pyproj.datadir.get_data_dir())
os.environ["PROJ_NETWORK"] = "OFF"
os.environ["PROJ_LIB"] = pyproj.datadir.get_data_dir()
os.environ["PROJ_DATA"] = pyproj.datadir.get_data_dir()

import rasterio

import sys
import logging

import dask
from dask.diagnostics import ProgressBar
import geopandas as gpd
import numpy as np
import planetary_computer
import pystac_client
import requests
import rioxarray  # noqa: F401 — registers .rio accessor on xarray objects
import xarray as xr
from odc.stac import load as odc_load

#import ndvi.python.NDVIDownloadPipeline.config as config
import config

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")
log = logging.getLogger(__name__)


# ── Phase 2 helpers ──────────────────────────────────────────────────────────

def load_aoi(shapefile_path: str):
    """Load shapefile and return GeoDataFrame in WGS84."""
    gdf = gpd.read_file(shapefile_path)
    if gdf.crs and gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    return gdf


def get_bbox(gdf: gpd.GeoDataFrame):
    """Return [west, south, east, north] from a GeoDataFrame."""
    return list(gdf.total_bounds)


def _configure_cdse_s3():
    """Configure GDAL/rasterio to read from CDSE's S3 endpoint."""
    config.ensure_cdse_credentials()
    os.environ["AWS_S3_ENDPOINT"] = config.CDSE_S3_ENDPOINT
    os.environ["AWS_ACCESS_KEY_ID"] = config.CDSE_S3_ACCESS_KEY
    os.environ["AWS_SECRET_ACCESS_KEY"] = config.CDSE_S3_SECRET_KEY
    os.environ["AWS_VIRTUAL_HOSTING"] = "FALSE"
    os.environ["AWS_HTTPS"] = "YES"
    # CDSE rate-limits aggressively — configure GDAL to retry on 429
    os.environ["GDAL_HTTP_MAX_RETRY"] = "5"
    os.environ["GDAL_HTTP_RETRY_DELAY"] = "2"
    log.info(f"GDAL S3 configured for CDSE endpoint: {config.CDSE_S3_ENDPOINT}")


def build_stac_client():
    """Return a pystac_client.Client configured for the selected backend."""
    if config.DATA_SOURCE == "mpc":
        return pystac_client.Client.open(
            config.MPC_STAC_URL,
            modifier=planetary_computer.sign_inplace,
        )
    elif config.DATA_SOURCE == "cdse":
        return pystac_client.Client.open(config.CDSE_STAC_URL)
    elif config.DATA_SOURCE == "aws":
        return pystac_client.Client.open(config.AWS_STAC_URL)
    else:
        raise ValueError(f"Unknown DATA_SOURCE: {config.DATA_SOURCE!r}")


def get_collection_name() -> str:
    collections = {
        "mpc": config.MPC_COLLECTION,
        "cdse": config.CDSE_COLLECTION,
        "aws": config.AWS_COLLECTION,
    }
    return collections[config.DATA_SOURCE]


_BAND_NAME_MAPPINGS = {
    "cdse": {"B04_10m": "B04", "B08_10m": "B08", "SCL_20m": "SCL"},
    "aws": {"red": "B04", "nir": "B08", "scl": "SCL"},
}


_PROJ_EXT_V1 = "https://stac-extensions.github.io/projection/v1.1.0/schema.json"
_PROJ_EXT_V2 = "https://stac-extensions.github.io/projection/v2.0.0/schema.json"


def _fixup_cdse_items(items):
    """Fix CDSE-specific STAC compatibility issues.

    1. Downgrade proj extension URI from v2 to v1 (odc-stac only supports v1).
    2. Translate proj:code → proj:epsg on each asset (v2 → v1 field name).
    3. Configure GDAL S3 for CDSE's endpoint (assets use s3://eodata/... HREFs).
    """
    _configure_cdse_s3()
    for item in items:
        # Fix extension URI
        if _PROJ_EXT_V2 in item.stac_extensions:
            item.stac_extensions.remove(_PROJ_EXT_V2)
        if _PROJ_EXT_V1 not in item.stac_extensions:
            item.stac_extensions.append(_PROJ_EXT_V1)

        for asset in item.assets.values():
            ef = asset.extra_fields

            # proj:code → proj:epsg
            if "proj:code" in ef and "proj:epsg" not in ef:
                code = ef.pop("proj:code")  # e.g. "EPSG:32630"
                if code is not None:
                    try:
                        ef["proj:epsg"] = int(code.split(":")[1])
                    except (IndexError, ValueError):
                        pass

    return items


def _ensure_proj_extension(items):
    """Inject the projection extension URI if missing.

    Some backends omit the proj extension declaration from their STAC
    items even though the data contains projection metadata.
    odc-stac requires it to be listed explicitly.
    """
    for item in items:
        if _PROJ_EXT_V1 not in item.stac_extensions:
            item.stac_extensions.append(_PROJ_EXT_V1)
    return items


def normalise_band_names(items):
    """Ensure asset keys match the canonical names in config.BANDS.

    Different backends use different asset key conventions. This renames
    them in-place so odc-stac can find the requested bands.
    """
    mapping = _BAND_NAME_MAPPINGS.get(config.DATA_SOURCE, {})
    if not mapping:
        return items
    for item in items:
        for old_key, new_key in mapping.items():
            if old_key in item.assets and new_key not in item.assets:
                item.assets[new_key] = item.assets.pop(old_key)
    return items


def search_catalog(client, bbox):
    """Query the STAC catalog and return a list of items."""
    query_filter = {}
    if config.MAX_CLOUD_COVER is not None:
        query_filter["eo:cloud_cover"] = {"lt": config.MAX_CLOUD_COVER}

    search = client.search(
        collections=[get_collection_name()],
        bbox=bbox,
        datetime=config.DATE_RANGE,
        query=query_filter if query_filter else None,
    )
    items = list(search.items())

    if not items:
        log.warning("STAC search returned 0 items — check date range, bbox, and cloud threshold.")
        sys.exit(0)

    log.info(f"STAC search returned {len(items)} items.")
    return items


# ── Phase 3 ──────────────────────────────────────────────────────────────────

def load_data(items, bbox) -> xr.Dataset:
    """Lazily load Sentinel-2 bands into an xarray Dataset via odc-stac."""
    if config.DATA_SOURCE == "cdse":
        items = _fixup_cdse_items(items)
    items = _ensure_proj_extension(items)
    if config.DATA_SOURCE in _BAND_NAME_MAPPINGS:
        items = normalise_band_names(items)

    ds = odc_load(
        items,
        bbox=bbox,
        bands=config.BANDS,
        resolution=config.RESOLUTION,
        groupby="solar_day",
        chunks=config.CHUNKS,
        dtype="float32",
        nodata=0,
        crs=config.OUTPUT_CRS,
        resampling={"SCL": "nearest", "*": "average"},
    )
    log.info(f"Loaded dataset: {ds.sizes}")
    return ds


def apply_cloud_mask(ds: xr.Dataset) -> xr.Dataset:
    """Mask pixels where SCL is not in VALID_SCL_CLASSES."""
    scl = ds["SCL"]
    valid = xr.zeros_like(scl, dtype=bool)
    for cls in config.VALID_SCL_CLASSES:
        valid = valid | (scl == cls)
    ds["B04"] = ds["B04"].where(valid)
    ds["B08"] = ds["B08"].where(valid)
    log.info("Cloud mask applied.")
    return ds


# ── Phase 4 ──────────────────────────────────────────────────────────────────

def compute_ndvi(ds: xr.Dataset) -> xr.DataArray:
    """Compute NDVI with division-by-zero guard and range clipping."""
    red = ds["B04"]
    nir = ds["B08"]
    denominator = nir + red
    ndvi = xr.where(denominator == 0, np.nan, (nir - red) / denominator)
    ndvi = ndvi.clip(-1, 1)
    ndvi.name = "NDVI"
    ndvi.attrs.update(
        units="dimensionless",
        long_name="NDVI",
        source="Sentinel-2 L2A",
        cloud_mask_classes=str(config.VALID_SCL_CLASSES),
    )
    # Preserve CRS from the source dataset (xr.where drops spatial attrs)
    ndvi = ndvi.rio.set_spatial_dims(x_dim=ds.rio.x_dim, y_dim=ds.rio.y_dim)
    ndvi.rio.write_crs(ds.rio.crs, inplace=True)
    
    log.info("NDVI computed.")
    return ndvi


def export(ndvi: xr.DataArray):
    """Export NDVI to the configured format(s)."""
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    formats = config.OUTPUT_FORMATS
    if isinstance(formats, str):
        formats = [formats]

    # Materialise the lazy Dask array with a progress bar (this is where
    # all the downloading and computation happens).
    log.info("Computing NDVI data cube (downloading pixels)...")
    with ProgressBar():
        ndvi = ndvi.compute()
    log.info("Compute complete.")

    for fmt in formats:
        if fmt == "netcdf":
            path = os.path.join(config.OUTPUT_DIR, "ndvi.nc")
            ndvi.to_netcdf(path)
            log.info(f"Exported NetCDF → {path}")

        elif fmt == "zarr":
            path = os.path.join(config.OUTPUT_DIR, "ndvi.zarr")
            ndvi.to_dataset(name="NDVI").to_zarr(path, mode="w")
            log.info(f"Exported Zarr → {path}")

        elif fmt == "cog":
            cog_dir = os.path.join(config.OUTPUT_DIR, "cog")
            os.makedirs(cog_dir, exist_ok=True)
            for t in ndvi.time.values:
                date_str = str(np.datetime_as_string(t, unit="D"))
                path = os.path.join(cog_dir, f"NDVI_{date_str}.tif")
                slice_ = ndvi.sel(time=t)
                slice_.rio.to_raster(path, driver="COG")
            log.info(f"Exported COGs → {cog_dir}/ ({len(ndvi.time)} files)")

        else:
            log.warning(f"Unknown output format: {fmt!r} — skipping.")


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Sentinel-2 NDVI data cube pipeline.",
    )
    parser.add_argument(
        "--source", choices=["mpc", "cdse", "aws"], default=None,
        help="Data source backend (default: from config.py)",
    )
    parser.add_argument(
        "--dates", default=None, metavar="START/END",
        help="Date range, e.g. 2023-01-01/2023-12-31 (default: from config.py)",
    )
    parser.add_argument(
        "--cloud-cover", type=float, default=None, metavar="PCT",
        help="Max scene-level cloud cover %%. Use -1 to disable (default: from config.py)",
    )
    parser.add_argument(
        "--no-cloud-mask", action="store_true",
        help="Disable pixel-level SCL cloud masking",
    )
    parser.add_argument(
        "--format", dest="formats", nargs="+",
        choices=["netcdf", "zarr", "cog"], default=None,
        help="Output format(s) (default: from config.py)",
    )
    parser.add_argument(
        "--output-dir", default=None, metavar="DIR",
        help="Output directory (default: from config.py)",
    )
    parser.add_argument(
        "--shape-path", default=None, metavar="FILE",
        help="shapefile directory (default: from config.py)",
    )
    return parser.parse_args()


def apply_cli_overrides(args):
    """Override config values with any CLI arguments provided."""
    if args.source is not None:
        config.DATA_SOURCE = args.source
    if args.dates is not None:
        config.DATE_RANGE = args.dates
    if args.cloud_cover is not None:
        config.MAX_CLOUD_COVER = None if args.cloud_cover < 0 else args.cloud_cover
    if args.no_cloud_mask:
        config.APPLY_CLOUD_MASK = False
    if args.formats is not None:
        config.OUTPUT_FORMATS = args.formats
    if args.output_dir is not None:
        config.OUTPUT_DIR = args.output_dir
    if args.shape_path is not None:
        config.SHAPEFILE_PATH = args.shape_path


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    apply_cli_overrides(args)

    # Set Dask concurrency (CDSE needs throttling due to rate limits)
    if config.DATA_SOURCE == "cdse":
        dask.config.set(num_workers=4)
    elif config.DASK_WORKERS is not None:
        dask.config.set(num_workers=config.DASK_WORKERS)
    log.info(f"Data source: {config.DATA_SOURCE}")
    log.info(f"Date range:  {config.DATE_RANGE}")

    # Phase 2
    gdf = load_aoi(config.SHAPEFILE_PATH)
    bbox = get_bbox(gdf)
    log.info(f"AOI bbox (WGS84): {bbox}")

    client = build_stac_client()
    items = search_catalog(client, bbox)

    # Phase 3
    ds = load_data(items, bbox)
    if config.APPLY_CLOUD_MASK:
        ds = apply_cloud_mask(ds)
    else:
        log.info("Cloud masking disabled — using raw reflectance.")

    # Phase 4
    ndvi = compute_ndvi(ds)
    export(ndvi)
    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()
