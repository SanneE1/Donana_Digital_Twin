# Sentinel-2 NDVI Data Cube Pipeline

Automated pipeline for generating cloud-masked NDVI data cubes from Sentinel-2 L2A imagery over any Area of Interest. Ships with a bundled example for the Doñana wetland, Spain.

## What it does

1. Reads an AOI from a shapefile (Doñana example included)
2. Searches a STAC catalog (MPC, CDSE, or AWS Earth Search) for Sentinel-2 L2A scenes
3. Loads only the required bands (Red, NIR, SCL) cropped to the AOI bounding box
4. Applies pixel-level cloud masking using the Scene Classification Layer (SCL)
5. Computes NDVI with division-by-zero protection and range clipping
6. Exports to NetCDF, Zarr, and/or Cloud-Optimized GeoTIFF
7. Generates validation plots (spatial maps, time series, histograms)

## Quick start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run with the bundled Doñana example

```bash
python pipeline.py
```

This uses the default config: Doñana AOI, AWS Earth Search, June 2023, NetCDF output. No accounts or credentials needed.

### 3. Use your own AOI

Edit `config.py`:

```python
SHAPEFILE_PATH = "/path/to/your/shapefile.shp"
OUTPUT_CRS = "EPSG:XXXXX"  # UTM zone covering your AOI
DATE_RANGE = "2023-01-01/2023-12-31"
```

### 4. Configure

| Parameter | Default | Description |
|---|---|---|
| `SHAPEFILE_PATH` | `examples/donana/AIS_Doñana.shp` | Path to AOI shapefile |
| `DATA_SOURCE` | `"aws"` | `"aws"` (Earth Search), `"mpc"` (Planetary Computer), or `"cdse"` (Copernicus) |
| `DATE_RANGE` | `"2023-06-01/2023-06-30"` | Temporal range in `YYYY-MM-DD/YYYY-MM-DD` format |
| `MAX_CLOUD_COVER` | `20` | Scene-level cloud cover threshold (%). Set to `None` to disable |
| `APPLY_CLOUD_MASK` | `True` | Pixel-level SCL masking. Set to `False` to use raw reflectance |
| `OUTPUT_FORMATS` | `["netcdf"]` | One or more of: `"netcdf"`, `"zarr"`, `"cog"` |
| `OUTPUT_CRS` | `"EPSG:32629"` | Output CRS (UTM 29N for Doñana; change for other AOIs) |
| `VALID_SCL_CLASSES` | `[2, 4, 5, 6, 7]` | SCL classes to keep (see below) |

### 5. CLI overrides

```bash
python pipeline.py --source aws             # override data source
python pipeline.py --dates 2022-01-01/2022-12-31  # override date range
python pipeline.py --cloud-cover 50         # allow cloudier scenes
python pipeline.py --cloud-cover -1         # disable scene-level cloud filter
python pipeline.py --no-cloud-mask          # disable pixel-level SCL masking
python pipeline.py --format netcdf zarr     # export to multiple formats
python pipeline.py --output-dir /tmp/ndvi   # custom output directory

python pipeline.py --shape_path /data/area.shp  # Added by me, hopefully it works
```

All arguments are optional — when omitted, values from `config.py` are used. Run `python pipeline.py --help` for the full list.

### 6. Visualise the results

```bash
python visualise.py               # load from NetCDF (default)
python visualise.py --from zarr   # load from Zarr store
```

Plots are saved to `output/plots/`.

## Cloud filtering

Cloud filtering operates at two levels, both independently configurable:

### Scene-level pre-filter (`MAX_CLOUD_COVER`)

The STAC query discards scenes where the overall cloud cover percentage exceeds this threshold. This reduces the number of scenes downloaded and speeds up processing.

- Set to `20` (default): only scenes with <20% cloud cover are fetched. This is standard in Sentinel-2 literature — it balances temporal coverage against download volume, since pixel-level masking handles the remaining clouds anyway.
- Increase to `50` or higher for cloudy regions or seasons where 20% discards too many scenes (e.g., Atlantic-facing coasts in winter).
- Set to `None`: all scenes are fetched regardless of cloud cover. Maximises temporal coverage at the cost of slower runs (more data to download). Pixel-level SCL masking still removes cloud pixels.

### Pixel-level mask (`APPLY_CLOUD_MASK`)

After loading, the SCL (Scene Classification Layer) band classifies each pixel. Only pixels in the `VALID_SCL_CLASSES` are kept; everything else becomes NaN.

| SCL Class | Label | Default |
|---|---|---|
| 2 | Dark Area Pixels | **Kept** — useful for wetlands (dark shallow water, saturated soil) |
| 4 | Vegetation | **Kept** |
| 5 | Bare Soils | **Kept** |
| 6 | Water | **Kept** — important for coastal/wetland AOIs |
| 7 | Unclassified | **Kept** — wet mudflats often fall here |
| 0, 1, 3, 8–11 | No data, defective, cloud shadow, clouds, etc. | **Masked → NaN** |

Adjust `VALID_SCL_CLASSES` in `config.py` for your study area. For example, if your AOI has no water bodies, you may want to remove class 6.

Set `APPLY_CLOUD_MASK = False` to disable pixel masking entirely (useful for inspecting raw data or debugging SCL behaviour).

## Spatial coverage note

When an AOI spans multiple Sentinel-2 MGRS tiles, some dates may have **partial spatial coverage** — one tile passes the cloud filter while another does not. This is normal for multi-tile analysis.

Implications:
- Each pixel may have a different number of valid observations across the time series
- The spatial mean NDVI per date may be biased toward whichever tiles were cloud-free
- Disabling the scene-level cloud filter (`MAX_CLOUD_COVER = None`) maximises spatial coverage at the cost of including cloudier scenes (pixel-level masking still removes cloud pixels)

## Data sources

### Rate limits & performance comparison

| | AWS Earth Search | MPC | CDSE (free tier) |
|---|---|---|---|
| **Auth for search** | None | None | Yes |
| **Auth for data** | None | SAS token (auto) | S3 keys |
| **Concurrent downloads** | High (S3) | High | **2** |
| **Bandwidth cap** | None | Not published | **~20 Mbps** |
| **Monthly volume** | None | Not published | **~1 TB** |
| **Data format** | COG | COG | JPEG2000 |
| **Practical speed** | Fast | **Fastest** (in our tests) | Slow |

CDSE is significantly slower than the other two backends due to strict rate limits and JPEG2000 format (no efficient range requests). The pipeline automatically throttles Dask and configures GDAL retries when using CDSE, but expect runs to take considerably longer.

### AWS Earth Search (Element 84) — default

- **No authentication required** — zero setup, works out of the box
- Full Sentinel-2 L2A archive (2015–present) as COGs on AWS
- Set `DATA_SOURCE = "aws"` in config (this is the default)
- Band names are normalised automatically (`red` → `B04`, `nir` → `B08`, `scl` → `SCL`)
- Ideal for new users, CI pipelines, or environments where installing `planetary-computer` is inconvenient

### Microsoft Planetary Computer (MPC) — fastest

- No account required for searching; SAS tokens are generated automatically via the `planetary-computer` package
- Full Sentinel-2 L2A archive (2015–present) as Cloud-Optimized GeoTIFFs
- Set `DATA_SOURCE = "mpc"` in config
- ~2x faster than AWS in our tests (Azure infrastructure may have better peering depending on your location)

### Copernicus Data Space Ecosystem (CDSE) — slow

- Requires a free account at [dataspace.copernicus.eu](https://dataspace.copernicus.eu)
- **Requires S3 access keys** (not OAuth2). Generate them at: [eodata-s3keysmanager.dataspace.copernicus.eu](https://eodata-s3keysmanager.dataspace.copernicus.eu/panel/s3-credentials)
- Set `DATA_SOURCE = "cdse"` in config and run the pipeline
- **Credentials are handled automatically.** On first run, the pipeline will:
  1. Check for `CDSE_S3_ACCESS_KEY` / `CDSE_S3_SECRET_KEY` environment variables
  2. If not found, check for a local credentials file at `.credentials/cdse.json`
  3. If neither exists, interactively prompt for your S3 keys, then save them to `.credentials/cdse.json` (permissions `600`, excluded from git)
- Band names are normalised automatically (e.g., `B04_10m` → `B04`)
- Data is read directly from CDSE's S3 endpoint (`eodata.dataspace.copernicus.eu`)
- **Expect slow performance** — free tier limits concurrency to 2 connections at ~20 Mbps, and data is stored as JPEG2000 (not COG)

You can also set credentials manually via environment variables if you prefer:
```bash
export CDSE_S3_ACCESS_KEY="your-access-key"
export CDSE_S3_SECRET_KEY="your-secret-key"
```

## Output formats

| Format | File(s) | Best for |
|---|---|---|
| **NetCDF** (default) | `output/ndvi.nc` | Single-file archival, QGIS, xarray reload |
| **Zarr** | `output/ndvi.zarr/` | Large cubes, parallel write, cloud storage |
| **COG** | `output/cog/NDVI_YYYY-MM-DD.tif` | Per-date sharing, tile servers, universal GIS |

Configure via `config.py` or CLI:
```bash
python pipeline.py --format netcdf              # default
python pipeline.py --format zarr                # Zarr only
python pipeline.py --format netcdf zarr cog     # all three in one run
```

### Reloading outputs in Python

```python
import xarray as xr

# NetCDF
ndvi = xr.open_dataarray("output/ndvi.nc")

# Zarr
ndvi = xr.open_zarr("output/ndvi.zarr")["NDVI"]

# COG (single date)
import rioxarray
ndvi = xr.open_dataarray("output/cog/NDVI_2023-06-15.tif", engine="rasterio")
```

### Notes

- **NetCDF** is a single self-contained file. Good default for most users.
- **Zarr** is a directory of chunks. Faster to write for large time series (Dask writes chunks in parallel) and works well with cloud storage (S3/GCS). Requires the `zarr` package.
- **COG** produces one georeferenced GeoTIFF per date with internal tiling and overviews. Opens in any GIS tool and can be served over HTTP. The CRS is embedded in each file.

## Visualisation outputs

| Plot | File | Purpose |
|---|---|---|
| Spot-check map | `output/plots/spot_check.png` | Verify spatial extent, CRS alignment, and cloud masking on a single date |
| Time series | `output/plots/time_series.png` | Spatial-mean NDVI with 10th–90th percentile band over time |
| Multi-panel | `output/plots/multi_panel.png` | Small-multiples grid of NDVI maps at evenly spaced dates |
| Histogram | `output/plots/histogram.png` | Per-date NDVI distribution heatmap — flags anomalous values |
| Interactive | `output/plots/interactive.html` | Slider-based map (requires `hvplot`) |

## Project structure

```
.
├── config.py              # All tuneable parameters
├── pipeline.py            # Main pipeline (AOI → STAC → load → mask → NDVI → export)
├── visualise.py           # Validation plots and data cube exploration
├── requirements.txt       # Python dependencies
├── examples/              # Bundled example AOIs (tracked in git)
│   └── donana/
│       └── AIS_Doñana.shp
├── .credentials/          # Local credential storage (not tracked in git)
│   └── cdse.json
├── data/                  # User data (not tracked in git)
└── output/                # Pipeline outputs (not tracked in git)
    ├── ndvi.nc
    ├── ndvi.zarr/
    ├── cog/
    └── plots/
```
