import os
from os import path
from os import getcwd
from os import makedirs
import sys

import tempfile
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import beta

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

import pyabc
from pyabc.transition import DiscreteJumpTransition, AggregatedTransition, MultivariateNormalTransition

try:
    # This works when running as a script
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
except NameError:
    # This works in Spyder/interactive environments where __file__ isn't defined
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(os.getcwd())))

sys.path.insert(0, parent_dir)

from python.ndvi_to_K_map_beta import ndvi_to_K_map

pyabc.settings.set_figure_params('pyabc')  # for beautified plots

# Start ------------------------------------------------------------------------------------------------

# prepare observation files ------------------------------------------------------------------------------------------------

# Donana 
obs_df = pd.read_csv(path.join(getcwd(), "data" , "original_data", "Rabbit_donana_KAI_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"), skipinitialspace=True)

obs_df = obs_df.melt(id_vars=['Fecha'],  # columns to keep as-is
                  value_vars=['Coto del Rey', 'Algaida-Sotos', 'Sabinar-Mogea', 
                              'RBD-este', 'Puntal', 'Marismillas', 'Abalario',
                              'Hinojos'],  # columns to pivot
                  var_name='transect',  # name for the new column with old column names
                  value_name='KAI') 

# get R functions
r = robjects.r
r['source'](path.join(getcwd(), 'R', 'function_rabbit_distances.R'))

get_sim_data = robjects.globalenv['get_sim_data']
csvToRaster = robjects.globalenv['csvToRaster']

hab_file = path.join(getcwd(), "data", "GIS_maps", "Rabbit_HabitatMap_500_Donana_Fordham_2013.asc")
transectsD = path.join(getcwd(), "data", 'original_data', "Rabbit_donana_KAI_PacoCarro", "Transect_oryctolagus.kml")
ndvi_map = path.join(getcwd(), "data", "model_input", "maps", "ndvi_monthly_maps")
flood_map = path.join(getcwd(), "data", "model_input", "maps", "LAST_floodmaps")

# parameters = {
#     "alpha": 2,
#     "beta_param": 2,
#     "K": 15,
#     "dens_opt": 5,
#     "R_lambda": 0.1,
#     "R_sigma": 0.1,
#     "obs_p": 0,
#     "obs_p_ndvi": 0}

# tmpdir = "temp_abc_folder"

# define model ------------------------------------------------------------------------------------------------
def model(parameters, n_repeats=3): 
    try:
        all_reps = []
        
        with tempfile.TemporaryDirectory() as tmpdir:
            hab_folder = path.join(tmpdir, "habitat_folder")
            makedirs(hab_folder, exist_ok=True)
        
            ndvi_to_K_map(ndvi_dir = ndvi_map,
                          flood_dir = flood_map,
                          output_map = hab_folder, 
                          alpha = parameters['alpha'], 
                          beta_param = parameters['beta_param'], 
                          K = parameters['K'],
                          start_year= 2003)
                    
            for r in range(n_repeats):
                sim_dir = path.join(tmpdir, f"sim_{r}")
                makedirs(sim_dir, exist_ok=True)
                
                cmd = [#"./Rabbit_solo_model",
                       path.join("Pascal_program","programs","Rabbit_solo_model.exe"), 
                       path.join("data", "model_input", "Parameters_historic_calibration.txt"), 
                       sim_dir,
                       str(parameters['dens_opt']), 
                       str(parameters['R_lambda']), 
                       str(parameters['R_sigma']),
                       path.join(getcwd(), "data", "model_input", "maps", "Rabbit_HabitatMap_500_Donana_Fordham_2013.txt"),
                       hab_folder
                       ] 
                mod = subprocess.run(cmd, capture_output=True, text=True, timeout=600) 
                if mod.returncode != 0: 
                    raise RuntimeError(f"Subprocess failed: {mod.stderr}")   
                
                obs_p = float(parameters['obs_p'])
                obs_p_ndvi = float(parameters['obs_p_ndvi'])
                
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    sim_data = get_sim_data(sim_dir, hab_file, transectsD, obs_df,
                                            obs_p, obs_p_ndvi)
                    
                    with (robjects.default_converter + pandas2ri.converter).context():
                        sim_data = robjects.conversion.get_conversion().rpy2py(sim_data)
                
                all_reps.append(sim_data)
                
            all_reps = pd.concat(all_reps, ignore_index=True)

        mean_sim_data = (
            all_reps
            .groupby(["transect", "Fecha"], as_index=False)
            .agg({"KAI": "mean"})
        )
        
        return {"sim_data": mean_sim_data}
    
    except Exception as e:
        print(f"Model run failed for parameters {parameters}: {e}")
        # re-raise so pyABC treats this as invalid and retries
        raise

# calculate model fit ------------------------------------------------------------------------------------------------
def distance(sim_data, _):
        
    obs_data = obs_df.copy()
    sim_data = sim_data['sim_data'].copy()
    date_col = 'Fecha'

    obs_data[date_col] = pd.to_datetime(obs_data[date_col], format='%d/%m/%Y')
    sim_data[date_col] = pd.to_datetime(sim_data[date_col], origin='1970-01-01', unit='D')
    
    merged = pd.merge(obs_data, sim_data, on=[date_col, "transect"], suffixes=('_obs', '_sim'), how='left')
    merged = merged.sort_values(date_col).reset_index(drop=True)
    
    se_per_obs = (merged['KAI_obs'] - merged['KAI_sim'])**2
    
    return(se_per_obs.sum())
        
# define priors ------------------------------------------------------------------------------------------------       
domain_K= np.arange(20, 50)
domain_densopt = np.arange(1,15)

transition = AggregatedTransition(
    mapping = {
        "alpha": MultivariateNormalTransition(),
        "beta_param": MultivariateNormalTransition(),
        "K": DiscreteJumpTransition(domain=domain_K, p_stay=0.5),
        "dens_opt": DiscreteJumpTransition(domain=domain_densopt, p_stay=0.5),
        "R_lambda": MultivariateNormalTransition(),
        "R_sigma": MultivariateNormalTransition(),
        'obs_p': MultivariateNormalTransition(),
        'obs_p_ndvi': MultivariateNormalTransition()
    })



prior = pyabc.Distribution(
    alpha = pyabc.RV("uniform", 1.1, 10),
    beta_param = pyabc.RV("uniform", 1.1, 10),
    K = pyabc.RV("rv_discrete", values=(domain_K, np.ones_like(domain_K) / len(domain_K))),
    dens_opt = pyabc.RV("rv_discrete", values=(domain_densopt, np.ones_like(domain_densopt) / len(domain_densopt))), 
    R_lambda = pyabc.RV("uniform", 0.01, 0.8),   
    R_sigma = pyabc.RV("uniform", 0.1, 5.0),
    obs_p = pyabc.RV("uniform", -4, 4),
    obs_p_ndvi = pyabc.RV("uniform", -4, 4)
    )            



# Run analysis ------------------------------------------------------------------------------------------------
sampler = pyabc.sampler.SingleCoreSampler()
abc = pyabc.ABCSMC(model, prior, distance, population_size=200, sampler = sampler,  
                   transitions=transition, eps = pyabc.QuantileEpsilon(alpha=0.5))    
        
# Initialize database for results
db_path = path.join("results", "rabbit_abc.db")
abc.new("sqlite:///" + db_path)

# Run ABC
history = abc.run(max_nr_populations=4)


# Save results ------------------------------------------------------------------------------------------------

# Save parameter estimates
df_final, w_final = history.get_distribution(m=0, t=history.max_t)

best_estimates = {}
all_params = ["alpha", "beta_param", "K",
              "dens_opt", "R_lambda", "R_sigma",
              'obs_p', 'obs_p_ndvi']

for param in all_params:
    if param in df_final.columns:
        # Calculate various statistics
        weighted_mean = np.average(df_final[param], weights=w_final)
        weighted_std = np.sqrt(np.average((df_final[param] - weighted_mean)**2, weights=w_final))
        
        best_estimates[param] = {
            'weighted_mean': weighted_mean,
            'weighted_std': weighted_std
            }

# Convert to DataFrame and save
best_estimates_df = pd.DataFrame.from_dict(best_estimates, orient='index')
best_estimates_df.index.name = 'parameter'
best_estimates_df.to_csv(path.join('results', 'parameter_estimates.csv'))



# Plot some of the results ------------------------------------------------------------------
fig, axes = plt.subplots(3, 3, figsize=(16, 12))
axes = axes.flatten()
all_param_ranges = {"alpha": (0, 5), 
                    "beta_param": (0, 5),
                    "K": (20, 50),
                    "dens_opt": (1, 15),
                    "R_lambda": (0.01, 0.8),
                    "R_sigma": (0.1, 5),
                    'obs_p': (-4, 4),
                    'obs_p_ndvi': (-4, 4)
                    }
                       
for i, param in enumerate(all_params):
    ax = axes[i]
    for t in range(history.max_t + 1):
        df, w = history.get_distribution(m=0, t=t)
        if param in df.columns:
            xmin, xmax = all_param_ranges[param]
            
            # Manual weighted KDE
            x_plot = np.linspace(xmin, xmax, 200)
            values = df[param].values
            
            # Create weighted KDE
            kde = stats.gaussian_kde(values, weights=w)
            y_plot = kde(x_plot)
            
            ax.plot(x_plot, y_plot, label=f"t={t}", alpha=0.7)
            
    ax.set_xlabel(param.replace('L_', ''))
    ax.set_title(param.replace('L_', ''))
    ax.set_xlim(all_param_ranges[param])
    if i == 0:
        ax.legend()

# Hide the extra subplot (9th one)
axes[8].axis('off')

plt.tight_layout()
plt.savefig(path.join('results', 'parameter_evolution.png'), dpi=300, bbox_inches='tight')
plt.show()

## ABC diagnostics
fig, arr_ax = plt.subplots(1, 3, figsize=(12, 4))

pyabc.visualization.plot_sample_numbers(history, ax=arr_ax[0])
pyabc.visualization.plot_epsilons(history, ax=arr_ax[1])
pyabc.visualization.plot_effective_sample_sizes(history, ax=arr_ax[2])

fig.tight_layout()
plt.savefig(path.join('results', 'parameter_abc_diagnostics.png'))


df_final, w_final = history.get_distribution(m=0, t=history.max_t)
print("\nFinal parameter estimates (weighted means):")
for param in all_params:
    if param in df_final.columns:
        weighted_mean = np.average(df_final[param], weights=w_final)
        print(f"{param}: {weighted_mean:.6f}")



# Save the final posterior distribution for forecasts -------------------------------------------------------------
df_final, w_final = history.get_distribution(m=0, t=history.max_t)

# Save particles and weights
df_final.to_csv(path.join('results', 'posterior_particles.csv'), index=False)
np.save(path.join('results', 'posterior_weights.npy'), w_final)

print(f"Saved {len(df_final)} posterior particles")



# NDVI relationship ------------------------------------------------------------------

indices = np.random.choice(len(df_final), size=500, replace=True, p=w_final)
sampled_params = df_final.iloc[indices].reset_index(drop=True)

# Create x values (NDVI range from 0 to 1)
x = np.linspace(0, 1, 1000)

# Show percentile envelope
plt.figure(figsize=(10, 6))

# Calculate beta PDF for all parameter combinations
all_curves = np.array([beta.pdf(x, row['alpha'], row['beta_param']) 
                       for _, row in df_final.iterrows()])

# Plot percentiles
plt.fill_between(x, np.percentile(all_curves, 2.5, axis=0), 
                 np.percentile(all_curves, 97.5, axis=0), 
                 alpha=0.3, label='95% CI')
plt.plot(x, np.mean(all_curves, axis=0), 'b--', linewidth=2, label='Mean')

plt.xlabel('NDVI', fontsize=12)
plt.ylabel('Suitability (Beta PDF)', fontsize=12)
plt.title('Beta Distribution Uncertainty', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(path.join('results', 'ndvi_beta_distribution.png'))
plt.show()


# Run historical simulations and plot -----------------------------------------------------

import analysis.run_history_sim

analysis.run_history_sim.main(N = 5)

# Plot results
robjects.r.assign("run_files", "results/history_summary")
robjects.r.source('R/plot_results.R')




