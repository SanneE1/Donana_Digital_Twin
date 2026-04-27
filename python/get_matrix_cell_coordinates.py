import rasterio
import numpy as np
from pyproj import Transformer
import csv
import os

def extract_and_transform_coordinates(raster_file, output_file, fromEPSG, toEPSG):
    """
    Extract coordinates from raster file and transform between EPSG
    
    Parameters:
    raster_file (str): Path to the raster file
    output_file (str): Path for the output coordinate file
    """
    
    # Create transformer from EPSG:3035 to EPSG:4326
    transformer = Transformer.from_crs(fromEPSG, toEPSG, always_xy=True)
    
    try:
        # Open the raster file
        with rasterio.open(raster_file) as src:
            print(f"Processing {raster_file}")
            print(f"Raster shape: {src.shape}")
            print(f"Raster CRS: {src.crs}")
            
            # Calculate total number of cells
            mat_dim = src.width * src.height
            print(f"Total cells to process: {mat_dim}")
            
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            # Open output file for writing
            with open(output_file, 'w', newline='') as outfile:
                writer = csv.writer(outfile, delimiter=' ')
                
                # Process cells in batches to manage memory
                batch_size = 10000
                processed = 0
                
                for start_idx in range(1, mat_dim + 1, batch_size):
                    end_idx = min(start_idx + batch_size - 1, mat_dim)
                    cell_indices = list(range(start_idx, end_idx + 1))
                    
                    # Convert cell indices to row, col
                    rows, cols = np.divmod(np.array(cell_indices) - 1, src.width)
                    
                    # Get coordinates using rasterio's transform
                    xs, ys = rasterio.transform.xy(src.transform, rows, cols)
                    
                    # Transform coordinates from EPSG:3035 to EPSG:4326
                    for x, y, row, col in zip(xs, ys, rows, cols):
                        try:
                            # Transform coordinates
                            lon, lat = transformer.transform(x, y)
                            
                            # Write longitude, latitude, col, row (matching R script output format)
                            writer.writerow([f"{lon:.6f}", f"{lat:.6f}", str(col + 1), str(row + 1)])
                            
                            processed += 1
                            
                            # Progress indicator
                            if processed % 50000 == 0:
                                print(f"Processed {processed}/{mat_dim} cells...")
                                
                        except Exception as e:
                            print(f"Warning: Could not transform coordinates for cell {start_idx + len(xs) - len(xs) + xs.index(x)}: {e}")
                            continue
                
                print(f"Completed processing {processed} cells")
                
    except Exception as e:
        print(f"Error processing file {raster_file}: {e}")
        raise

