import sys
import os
import glob

def transform_asc_file(input_path, output_path):
    """
    Transform an ASC file by replacing the standard 6-line header with a simplified
    single-line 'ncols nrows' header.
    
    Args:
        input_path (str): Path to the input ASC file
        output_path (str): Path for the output file. If None, will create a file with same name
    
    Returns:
        str: Path to the output file
    """
    
    ncols = None
    nrows = None
    
    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            # Read and parse the header
            for i in range(6):  # Standard ASC has 6 header lines
                line = infile.readline().strip()
                parts = line.split()
                
                if i == 0:  # ncols line
                    ncols = parts[1]
                elif i == 1:  # nrows line
                    nrows = parts[1]
            
            if ncols is None or nrows is None:
                raise ValueError("Could not find ncols or nrows in the ASC header")
            
            # Write the simplified header
            outfile.write(f"{ncols} {nrows}\n")
            
            # Copy the rest of the file (the raster data)
            for line in infile:
                outfile.write(line.rstrip() + '\n')
                
        print(f"Successfully transformed {input_path} to {output_path}")
        return output_path
        
    except Exception as e:
        print(f"Error transforming ASC file: {e}")
        if os.path.exists(output_path):
            try:
                os.remove(output_path)  # Clean up partial output file
            except:
                pass
        return None

def process_folder(asc_folder, output_folder):
    """
    Process all .asc files in the specified folder and subfolders, creating corresponding .txt files
    while maintaining the same directory structure.
    
    Args:
        asc_folder (str): Path to the folder containing .asc files
        output_folder (str): Path to where the .txt files need to be printed
    """
    # Make sure the source folder exists
    if not os.path.isdir(asc_folder):
        print(f"Error: Folder '{asc_folder}' does not exist")
        return
    
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    # Use os.walk to recursively find all .asc files
    asc_files = []
    for root, dirs, files in os.walk(asc_folder):
        for file in files:
            if file.lower().endswith('.asc'):
                asc_files.append(os.path.join(root, file))
    
    if not asc_files:
        print(f"No .asc files found in '{asc_folder}' or its subfolders")
        return
    
    print(f"Found {len(asc_files)} .asc files to process")
    
    # Process each file
    for asc_file in asc_files:
        # Calculate the relative path from the source folder
        rel_path = os.path.relpath(asc_file, asc_folder)
        
        # Get the directory structure and filename
        rel_dir = os.path.dirname(rel_path)
        filename = os.path.basename(rel_path)
        
        # Create the corresponding output directory structure
        output_dir = os.path.join(output_folder, rel_dir) if rel_dir else output_folder
        os.makedirs(output_dir, exist_ok=True)
        
        # Create output path with .txt extension
        base_name, _ = os.path.splitext(filename)
        output_filename = f"{base_name}.txt"
        output_file = os.path.join(output_dir, output_filename)
        
        # Transform the file
        transform_asc_file(asc_file, output_file)
