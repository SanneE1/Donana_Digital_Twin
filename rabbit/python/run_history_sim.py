# Forecast with calibrated parameters
import pandas as pd
import numpy as np
import os
import sys

from os import path
from os import getcwd
from os import makedirs

import tempfile
import subprocess

# Add the parent directory (Rabbit_forecast_DT) to the path (so it can be found when running in anaconda)
try:
    # This works when running as a script
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
except NameError:
    # This works in Spyder/interactive environments where __file__ isn't defined
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(os.getcwd())))
sys.path.insert(0, parent_dir)

from python.ndvi_to_K_map import ndvi_to_K_map
import python.get_summary_maps_from_simulations as summary_maps

def main(output_dir= path.join('results', 'history_summary'), 
         N=100, replicates=10):
    # Load info/files
    ndvi_map = path.join(getcwd(), "data", "model_input", "maps", "ndvi_monthly_maps")
    flood_map = path.join(getcwd(), "data", "model_input", "maps", "LAST_floodmaps")
    
    # Load the posterior
    df_posterior = pd.read_csv(path.join('results', 'posterior_particles.csv'))
    w_posterior = np.load(path.join('results', 'posterior_weights.npy'))
    
    # Draw N parameter sets from the posterior
    # N = 100  # number of forecast runs
    # replicates = 10
    indices = np.random.choice(len(df_posterior), size=N, replace=True, p=w_posterior)
    sampled_params = df_posterior.iloc[indices].reset_index(drop=True)
    
    # idx = 1
    # params = sampled_params.iloc[idx]
    # tmpdir = 'temp_folder'
    
    # Now run forecasts with each parameter set
    with tempfile.TemporaryDirectory() as tmpdir: 
        for idx, params in sampled_params.iterrows():
            print(f"parameter sample: {idx}")
            
            alpha = params['alpha']
            beta_param = params['beta_param']
            K = int(params['K'])
            dens_opt = int(params['dens_opt'])
            R_lambda = params['R_lambda']
            R_sigma = params['R_sigma']
            obs_p = params['obs_p']
            obs_p_ndvi = params['obs_p_ndvi']
            
            for rep in range(0, replicates):
                sim_folder = path.join(tmpdir, str(idx) + str(rep))
                hab_folder = path.join(sim_folder, "habitat_folder")
                makedirs(sim_folder, exist_ok=True)
                makedirs(hab_folder, exist_ok=True)
            
                ndvi_to_K_map(ndvi_dir = ndvi_map,
                              flood_dir = flood_map,
                              output_map = hab_folder, 
                              alpha = alpha, 
                              beta_param = beta_param, 
                              K = K,
                              start_year = 2010)
            
                cmd = [path.join("Pascal_program","programs","Rabbit_solo_model.exe"), 
                       path.join("data", "model_input", "Parameters_historic_simulation.txt"), 
                       sim_folder,
                       str(dens_opt), 
                       str(R_lambda), 
                       str(R_sigma),
                       path.join(hab_folder, '2010_01.txt'),
                       hab_folder
                       ] 
                mod = subprocess.run(cmd, capture_output=True, text=True, timeout=500)
        
        summary_maps.main(tmpdir, output_dir)
    


if __name__ == "__main__":
    main()

