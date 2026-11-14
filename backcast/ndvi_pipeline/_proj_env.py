"""Force PROJ to use rasterio's bundled data files.

Some Python environments end up with mismatched PROJ databases across the
geospatial stack (pyproj, rasterio, fiona each ship their own copy). Symptoms
range from spurious "database disk image is malformed" errors to silent CRS
lookup failures.

Setting PROJ_DATA / PROJ_LIB to rasterio's bundled directory resolves the
conflict. Existing env-var values are respected (setdefault), so a user who
has deliberately pointed PROJ elsewhere is not overridden.

Import this module BEFORE any other geospatial library (geopandas, rioxarray,
pyproj, ...) so PROJ picks up the env vars on first use.
"""
import os

try:
    import rasterio as _rasterio
except ImportError:
    pass
else:
    _proj_data = os.path.join(os.path.dirname(_rasterio.__file__), "proj_data")
    if os.path.isdir(_proj_data):
        os.environ.setdefault("PROJ_DATA", _proj_data)
        os.environ.setdefault("PROJ_LIB", _proj_data)
