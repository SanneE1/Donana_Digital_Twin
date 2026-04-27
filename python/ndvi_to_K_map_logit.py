import numpy as np
import math
import warnings
from pathlib import Path

def ndvi_to_K_map(ndvi_dir, flood_dir, output_map, alpha, beta_param, K, start_year=None):
    files = list(Path(ndvi_dir).glob('*'))
    flood_files = list(Path(flood_dir).glob('*'))
    
    # Create a dictionary mapping dates to flood file paths for quick lookup
    flood_dict = {}
    for flood_file in flood_files:
        # Extract date identifier from flood filename (e.g., 'YYYY_MM' from 'YYYY_MM.txt')
        date_key = flood_file.stem
        flood_dict[date_key] = flood_file
    
    # Filter files by minimum year if specified
    if start_year is not None:
        filtered_files = []
        for file in files:
            try:
                # Extract year from filename (assuming format YYYY_MM.txt)
                year = int(file.stem.split('_')[0])
                if year >= start_year:
                    filtered_files.append(file)
            except (ValueError, IndexError):
                # Skip files that don't match the expected format
                continue
        files = filtered_files 
        
    for file in files:
        # Read first line (header)
        with open(file, 'r') as f:
            first_line = f.readline()
        
        # Read data (skip first line)
        df = np.loadtxt(file, skiprows=1)
        
        # Check if flood data exists for this time period
        date_key = file.stem
        flood_mask = None
        if date_key in flood_dict:
            # Load flood data (skip header if it has one)
            flood_data = np.loadtxt(flood_dict[date_key], skiprows=1)
            # Create mask: True where NOT flooded (0), False where flooded (1)
            flood_mask = (flood_data == 0)
        
        # Handle invalid values (e.g., -9999 for no-data)
        # Mask invalid values before transformation
        valid_mask = np.isfinite(df) & (df >= 0)
        
        # Combine valid NDVI mask with flood mask if flood data exists
        if flood_mask is not None:
            valid_mask = valid_mask & flood_mask
        
        # Initialize output array with zeros or a no-data value
        df_K = np.zeros_like(df, dtype=int)
        
        # Only transform valid NDVI values
        if np.any(valid_mask):
            
            # Apply beta distribution only to valid values
            df_b = np.zeros_like(df)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                df_b[valid_mask] = 1/(1+math.exp(-(alpha + beta_param * df[valid_mask])))
            
            # Scale by K and round
            df_K_temp = np.zeros_like(df)
            df_K_temp[valid_mask] = df_b[valid_mask] * K
            
            # Round and convert to int only for valid values
            df_K[valid_mask] = np.round(df_K_temp[valid_mask]).astype(int)
            
            # Set invalid values to 0 or a specific no-data value
            df_K[~valid_mask] = 0
        
        # Write output
        output_file = Path(output_map) / file.name
        with open(output_file, 'w') as f:
            f.write(first_line)
            np.savetxt(f, df_K, fmt='%d', delimiter=' ')