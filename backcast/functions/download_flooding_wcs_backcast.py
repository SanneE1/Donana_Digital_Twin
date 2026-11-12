# Backcast version of environmental_data/functions/download_flooding_wcs.py: skips maps already on disk.
#
#   python download_flooding_wcs_backcast.py [output_folder]
from owslib.wcs import WebCoverageService
import os
import sys

# WCS endpoint
wcs_url = "https://geoserver.icts-donana.es/geoserver/inundacionv1/wcs?"

# Folder to save downloaded rasters
output_folder = sys.argv[1] if len(sys.argv) > 1 else "data/backcast/LAST_maps"
os.makedirs(output_folder, exist_ok=True)

# Connect to WCS
wcs = WebCoverageService(wcs_url, version='2.0.1')

# Loop over all coverages
for coverage_id in wcs.contents:
    out_path = os.path.join(output_folder, f"{coverage_id}.tif")
    if os.path.exists(out_path):
        continue

    print(f"Downloading {coverage_id}...")
    
    # Get the coverage as GeoTIFF
    response = wcs.getCoverage(
        identifier=coverage_id,
        format='image/tiff;application=geotiff'
    )
    
    # Save to file
    with open(out_path, 'wb') as f:
        f.write(response.read())
    
    print(f"Saved to {out_path}")

print("All coverages downloaded successfully!")

