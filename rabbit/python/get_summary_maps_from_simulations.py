import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import sys

# Use bottleneck if available for faster operations
try:
    import bottleneck as bn
    USE_BOTTLENECK = True
    print("Using Bottleneck for accelerated computations")
except ImportError:
    USE_BOTTLENECK = False
    print("Bottleneck not found, using NumPy (consider installing: pip install bottleneck)")

def calculate_matrix_statistics(folder_paths, output_folder, confidence_level=0.95):
    """
    Calculate mean and 95% CI for each cell across multiple CSV matrices.
    
    Parameters:
    -----------
    folder_paths : list of str
        List of folder paths containing the CSV matrices
    output_folder : str
        Path where summary files will be saved
    confidence_level : float
        Confidence level for the interval (default: 0.95)
    """
    
    # Create output folder if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get list of CSV files from the first folder (assuming all folders have same files)
    first_folder = Path(folder_paths[0], "maps")
    csv_files = sorted(first_folder.glob("*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in {first_folder}")
        return
    
    print(f"Found {len(csv_files)} CSV files")
    print(f"Processing {len(folder_paths)} folders")
    
    # Process each CSV file
    for csv_file in csv_files:
        filename = csv_file.name
        #print(f"\nProcessing {filename}...")
        
        # Collect all matrices for this file
        matrices = []
        
        for folder in folder_paths:
            file_path = Path(folder, "maps") / filename
            if file_path.exists():
                df = pd.read_csv(file_path, header=None)
                matrices.append(df.values)
            else:
                print(f"  Warning: {file_path} not found, skipping")
        
        if not matrices:
            print(f"  No valid matrices found for {filename}, skipping")
            continue
        
        # Stack matrices along a new axis (now shape is: n_folders x rows x cols)
        matrices_array = np.stack(matrices, axis=0)
        n_samples = matrices_array.shape[0]
        
        #print(f"  Computing statistics across {n_samples} matrices...")
        
        # Calculate mean across folders (axis=0)
        mean_matrix = np.mean(matrices_array, axis=0)
        
        # Calculate standard error
        std_matrix = np.std(matrices_array, axis=0, ddof=1)
        se_matrix = std_matrix / np.sqrt(n_samples)
        
        # Calculate 95% CI using t-distribution
        # For large n, this approaches the normal distribution
        t_value = stats.t.ppf((1 + confidence_level) / 2, n_samples - 1)
        margin_error = t_value * se_matrix
        
        lower_ci = mean_matrix - margin_error
        upper_ci = mean_matrix + margin_error
        
        # Save results
        base_name = filename.replace('.csv', '')
        
        pd.DataFrame(mean_matrix).to_csv(
            output_path / f"{base_name}_mean.csv", 
            index=False, 
            header=False
        )
        pd.DataFrame(std_matrix).to_csv(
            output_path / f"{base_name}_std.csv",
            index=False,
            header=False
        )
        pd.DataFrame(lower_ci).to_csv(
            output_path / f"{base_name}_lower_ci.csv", 
            index=False, 
            header=False
        )
        pd.DataFrame(upper_ci).to_csv(
            output_path / f"{base_name}_upper_ci.csv", 
            index=False, 
            header=False
        )
        


def main(parent_directory, output_folder):
    """
    Process all subdirectories and calculate matrix statistics.
    
    Parameters:
    -----------
    parent_directory : str or Path
        Parent directory containing subdirectories to process
    output_folder : str or Path
        Output folder for summary statistics
    """
    parent_dir = Path(parent_directory)
    output_folder = Path(output_folder)
    
    # Get all subdirectories in the parent directory
    folders = [str(f) for f in parent_dir.iterdir() if f.is_dir()]
    
    if not folders:
        print(f"No subdirectories found in {parent_dir}")
        return None  # or raise an error if you prefer
    
    print(f"Found {len(folders)} folders to process:")
    for folder in folders:
        print(f"  - {folder}")
    
    # Run the calculation
    result = calculate_matrix_statistics(folders, output_folder)
    return result


if __name__ == "__main__":
    # Usage: python script.py <parent_directory> <output_folder>
    # Example: python script.py ./data_folders ./summary_statistics

    parent_dir = sys.argv[1]
    output_folder = sys.argv[2]
    
    main(parent_dir, output_folder)
    
 