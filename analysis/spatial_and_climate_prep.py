import os
from python.transform_asc_to_input_txt_map import process_folder
from python.get_matrix_cell_coordinates import extract_and_transform_coordinates

from owslib.wcs import WebCoverageService

os.environ["R_HOME"] =  "C:\Program Files\R\R-4.5.1"  
import rpy2.robjects as robjects


# Get CORINE map ----------------------------------------------------------------------------------
robjects.r.assign("corine_raster_path", "data/original_data/U2018_CLC2018_V2020_20u1.tif")
robjects.r.assign("output_folder", "data/GIS_maps/")
robjects.r.assign("DT_borders_path", "data/original_data/Donana_DT_border/Limite_Don╠âana.shp")

# Get CORINE landcover type map
robjects.r.source('R/Format_landcover_map.R')

# Get historical climate data ----------------------------------------------------------------------------------

#firs central cell coordinates
rast_file_DT = "data/GIS_maps/Rabbit_HabitatMap_500_Donana_Fordham_2013.asc"
output_DT = os.path.join("data", "coordinates_DT_EPSG4326.txt")
extract_and_transform_coordinates(rast_file_DT, output_DT, "EPSG:25829", "EPSG:4326")

# Get historical climate ----------------------------------------------------------------------------------

robjects.r.source('R/Format_historical_climate_rabbit.R')


# Download EBD_LAST Flooding maps --------------------------------------------------------------------------------------------
# WCS endpoint
wcs_url = "https://geoserver.icts-donana.es/geoserver/inundacionv1/wcs?"

# Folder to save downloaded rasters
output_folder = "data/original_data/LAST_maps/Inundacion/"
os.makedirs(output_folder, exist_ok=True)

# Connect to WCS
wcs = WebCoverageService(wcs_url, version='2.0.1')

# Loop over all coverages
for coverage_id in wcs.contents:
    print(f"Downloading {coverage_id}...")
    
    # Get the coverage as GeoTIFF
    response = wcs.getCoverage(
        identifier=coverage_id,
        format='image/tiff;application=geotiff'
    )
    
    # Save to file
    out_path = os.path.join(output_folder, f"{coverage_id}.tif")
    with open(out_path, 'wb') as f:
        f.write(response.read())
    
    print(f"Saved to {out_path}")

print("All coverages downloaded successfully!")


# Create monthly maps for ndvi and flooding --------------------------------------------------------------------------------------------
robjects.r.source('R/Format_landcover_map')


# Get ASC maps to input format ----------------------------------------------------------------------------------
process_folder("data/GIS_maps/", "data/model_input/maps")
