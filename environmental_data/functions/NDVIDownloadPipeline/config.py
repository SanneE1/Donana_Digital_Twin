"""
Pipeline configuration for Sentinel-2 NDVI data cube generation.
"""
import json
import os

_PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── AOI ──────────────────────────────────────────────────────────────────────
# Path to the shapefile defining the Area of Interest.
# The bundled Doñana example is used by default.
_EXAMPLE_DIR = os.path.join(_PROJECT_DIR, "examples", "donana")
SHAPEFILE_PATH = next(
    os.path.join(_EXAMPLE_DIR, f)
    for f in os.listdir(_EXAMPLE_DIR)
    if f.endswith(".shp")
)

# ── Temporal range ───────────────────────────────────────────────────────────
DATE_RANGE = "2023-06-01/2023-06-30"

# ── Cloud cover threshold (%) for STAC pre-filter ───────────────────────────
# Set to None to disable scene-level cloud filtering (fetch all scenes).
MAX_CLOUD_COVER = 20

# ── Pixel-level cloud masking ───────────────────────────────────────────────
# When True, the SCL band is used to mask clouds/shadows per pixel.
# When False, raw reflectance is used as-is (useful for inspecting SCL behaviour).
APPLY_CLOUD_MASK = True

# ── Data source: "mpc", "cdse", or "aws" ─────────────────────────────────────
DATA_SOURCE = "aws"

# MPC settings
MPC_STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
MPC_COLLECTION = "sentinel-2-l2a"

# CDSE settings
CDSE_STAC_URL = "https://catalogue.dataspace.copernicus.eu/stac"
CDSE_COLLECTION = "sentinel-2-l2a"
CDSE_S3_ENDPOINT = "eodata.dataspace.copernicus.eu"

# AWS Earth Search settings
AWS_STAC_URL = "https://earth-search.aws.element84.com/v1"
AWS_COLLECTION = "sentinel-2-l2a"

_CREDENTIALS_DIR = os.path.join(_PROJECT_DIR, ".credentials")
_CDSE_CREDENTIALS_FILE = os.path.join(_CREDENTIALS_DIR, "cdse.json")


def _load_cdse_credentials():
    """Load CDSE S3 credentials from env vars, local file, or interactive prompt.

    CDSE requires S3 access keys (not OAuth2) for reading raster data.
    Generate them at: https://dataspace.copernicus.eu → Dashboard → S3 Access Keys.

    Lookup order:
      1. Environment variables CDSE_S3_ACCESS_KEY / CDSE_S3_SECRET_KEY
      2. Local credentials file (.credentials/cdse.json)
      3. Interactive prompt (saves to .credentials/cdse.json for next time)
    """
    # 1. Environment variables
    access_key = os.environ.get("CDSE_S3_ACCESS_KEY", "")
    secret_key = os.environ.get("CDSE_S3_SECRET_KEY", "")
    if access_key and secret_key:
        return access_key, secret_key

    # 2. Local credentials file
    if os.path.exists(_CDSE_CREDENTIALS_FILE):
        with open(_CDSE_CREDENTIALS_FILE) as f:
            creds = json.load(f)
        access_key = creds.get("s3_access_key", "")
        secret_key = creds.get("s3_secret_key", "")
        if access_key and secret_key:
            return access_key, secret_key

    # 3. Interactive prompt
    print("CDSE S3 credentials not found.")
    print("Generate S3 access keys at:")
    print("  https://eodata-s3keysmanager.dataspace.copernicus.eu/panel/s3-credentials")
    print()
    access_key = input("CDSE S3 Access Key: ").strip()
    secret_key = input("CDSE S3 Secret Key: ").strip()
    if not access_key or not secret_key:
        raise ValueError("CDSE S3 credentials are required when DATA_SOURCE='cdse'.")

    # Save for next time
    os.makedirs(_CREDENTIALS_DIR, exist_ok=True)
    with open(_CDSE_CREDENTIALS_FILE, "w") as f:
        json.dump({"s3_access_key": access_key, "s3_secret_key": secret_key}, f, indent=2)
    os.chmod(_CDSE_CREDENTIALS_FILE, 0o600)
    print(f"Credentials saved to {_CDSE_CREDENTIALS_FILE}")
    return access_key, secret_key


# CDSE credentials — loaded lazily so CLI --source overrides work.
CDSE_S3_ACCESS_KEY = ""
CDSE_S3_SECRET_KEY = ""


def ensure_cdse_credentials():
    """Load CDSE S3 credentials if not already set. Called by the pipeline when needed."""
    global CDSE_S3_ACCESS_KEY, CDSE_S3_SECRET_KEY
    if not CDSE_S3_ACCESS_KEY or not CDSE_S3_SECRET_KEY:
        CDSE_S3_ACCESS_KEY, CDSE_S3_SECRET_KEY = _load_cdse_credentials()

# ── Bands & masking ─────────────────────────────────────────────────────────
BANDS = ["B04", "B08", "SCL"]
VALID_SCL_CLASSES = [2, 4, 5, 6, 7]

# ── Loading parameters ──────────────────────────────────────────────────────
RESOLUTION = 10  # metres
OUTPUT_CRS = "EPSG:32629"  # UTM zone 29N (covers Doñana example; change for other AOIs)
CHUNKS = {"time": 1, "x": 2048, "y": 2048}

# ── Concurrency ───────────────────────────────────────────────────────────
# Number of Dask threads for downloading/computing.
# None = use Dask default (one thread per CPU core).
# CDSE is automatically capped at 4 regardless of this setting.
DASK_WORKERS = None

# ── Output ───────────────────────────────────────────────────────────────────
# One or more of: "netcdf", "zarr", "cog"
OUTPUT_FORMATS = ["netcdf"]
OUTPUT_DIR = os.path.join(_PROJECT_DIR, "output")
