import numpy as np
from matplotlib.patches import Patch
from settings import *
import os

def plot_mean_errors(smm_dir, ode_dir, temporal_hybrid_dir, spatial_hybrid_dir, save_dir, tmax, comp_index, num_regions, figsize):
    print("Mean errors:")
    errors_ode = []
    errors_temporal = []
    errors_spatial = []
    fig, ax = plt.subplots(figsize=figsize)
    model_colors = [colors["dark grey"], colors["middle blue"], colors["middle green"]]
    for r in range(num_regions):
        mean_smm = pd.read_csv(smm_dir + f"/means.csv")
        mean_smm = mean_smm[mean_smm.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)]
        # ODE
        mean_ode = pd.read_csv(ode_dir + f"/means.csv")
        mean_ode = mean_ode[mean_ode.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)]
        err_ode = np.square(abs(mean_smm - mean_ode)/1.)
        max_err_ode = np.max(abs(mean_smm - mean_ode)/1.)
        errors_ode.append(np.mean(err_ode))
        # Temporal-hybrid model
        mean_temporal_hybrid = pd.read_csv(temporal_hybrid_dir + f"/means.csv")
        mean_temporal_hybrid = mean_temporal_hybrid[mean_temporal_hybrid.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)]
        err_temp_hybrid = np.square(abs(mean_smm - mean_temporal_hybrid)/1.)
        max_err_temp_hybrid = np.max(abs(mean_smm - mean_temporal_hybrid)/1.)
        errors_temporal.append(np.mean(err_temp_hybrid))
        # Spatial-hybrid model
        mean_spatial_hybrid = pd.read_csv(spatial_hybrid_dir + f"/{r}_means.csv")
        mean_spatial_hybrid = mean_spatial_hybrid[mean_spatial_hybrid.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)] 
        err_spatial_hybrid = np.square(abs(mean_smm - mean_spatial_hybrid)/1.)
        max_err_spatial_hybrid = np.max(abs(mean_smm - mean_spatial_hybrid)/1.)
        errors_spatial.append(np.mean(err_spatial_hybrid))
        print("Region", r, "ODE:", "MSE ", errors_ode[-1], "Max AE ", max_err_ode, "Temporal-hybrid:", "MSE ", errors_temporal[-1], "Max AE ", max_err_temp_hybrid, "Spatial-hybrid:", "MSE ", errors_spatial[-1], "Max AE ", max_err_spatial_hybrid)
        
    x = np.arange(num_regions)
    bar_width = 0.25
    ax.bar(x - bar_width, errors_ode, color=model_colors[0], width=bar_width, label="MoM model")
    ax.bar(x, errors_temporal, color=model_colors[1], width=bar_width, label="Temporal-hybrid model")
    ax.bar(x + bar_width, errors_spatial, color=model_colors[2], width=bar_width, label="Spatial-temporal-hybrid model")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Region {i}" for i in range(num_regions)])
    ax.set_ylabel(r"MSPE($\mu_I$)")
    plt.tight_layout()
    fig.savefig(f"{save_dir}/mean_error_all_regions.png", dpi=dpi)
    
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/all_regions_legend.png", dpi=dpi, bbox_inches='tight')
        
def plot_var_errors(smm_dir, ode_dir, temporal_hybrid_dir, spatial_hybrid_dir, save_dir, tmax, comp_index, num_regions, figsize):
    print("Variance errors:")
    errors_ode = []
    errors_temporal = []
    errors_spatial = []
    fig, ax = plt.subplots(figsize=figsize)
    model_colors = [colors["dark grey"], colors["middle blue"], colors["middle green"]]
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        mean_smm = pd.read_csv(smm_dir + f"/moments.csv")
        mean_smm = mean_smm[mean_smm.Time <= tmax][col_name].iloc[:]
        # ODE
        mean_ode = pd.read_csv(ode_dir + f"/moments.csv")
        mean_ode = mean_ode[mean_ode.Time <= tmax][col_name].iloc[:]
        err_ode = np.square(abs(mean_smm - mean_ode)/1.)
        max_err_ode = np.max(abs(mean_smm - mean_ode)/1.)
        errors_ode.append(np.mean(err_ode))
        # Temporal-hybrid model
        mean_temporal_hybrid = pd.read_csv(temporal_hybrid_dir + f"/moments.csv")
        mean_temporal_hybrid = mean_temporal_hybrid[mean_temporal_hybrid.Time <= tmax][col_name].iloc[:]
        err_temp_hybrid = np.square(abs(mean_smm - mean_temporal_hybrid)/1.)
        max_err_temp_hybrid = np.max(abs(mean_smm - mean_temporal_hybrid)/1.)
        errors_temporal.append(np.mean(err_temp_hybrid))
        # Spatial-hybrid model
        mean_spatial_hybrid = pd.read_csv(spatial_hybrid_dir + f"/{r}_moments.csv")
        mean_spatial_hybrid = mean_spatial_hybrid[mean_spatial_hybrid.Time <= tmax][col_name].iloc[:]
        err_spatial_hybrid = np.square(abs(mean_smm - mean_spatial_hybrid)/1.)
        max_err_spatial_hybrid = np.max(abs(mean_smm - mean_spatial_hybrid)/1.)
        errors_spatial.append(np.mean(err_spatial_hybrid))
        print("Max variance in SMM:", np.max(mean_smm))
        print("Region", r, "ODE:", "MSE ", errors_ode[-1], "Max AE ", max_err_ode, "Temporal-hybrid:", "MSE ", errors_temporal[-1], "Max AE ", max_err_temp_hybrid, "Spatial-hybrid:", "MSE ", errors_spatial[-1], "Max AE ", max_err_spatial_hybrid)
    
    x = np.arange(num_regions)
    bar_width = 0.25
    ax.bar(x - bar_width, errors_ode, color=model_colors[0], width=bar_width, label="MoM model")
    ax.bar(x, errors_temporal, color=model_colors[1], width=bar_width, label="Temporal-hybrid model")
    ax.bar(x + bar_width, errors_spatial, color=model_colors[2], width=bar_width, label="Spatial-temporal-hybrid model")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Region {i}" for i in range(num_regions)])
    ax.set_ylabel(r"MSPE($\sigma_I^2$)")
    plt.tight_layout()
    fig.savefig(f"{save_dir}/var_error_all_regions.png", dpi=dpi)
    
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/all_regions_legend.png", dpi=dpi, bbox_inches='tight')
    
if __name__ == "__main__":
    figsize = (5, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel/Comparion"
    temporal_hybrid_model = "Hybrid2"
    spatial_temporal_hybrid_model = "Spatial-Hybrid2"
    ode_model = "Moments"
    config = "config_2r_10_100_k1"
    num_regions = 2
    condition = "var_gradient"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    tmin = 0
    tmax = 90
    
    smm_dir = f"{dir}/SMM/{config}"
    temporal_hybrid_dir = f"{dir}/{temporal_hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    spatial_hybrid_dir = f"{dir}/{spatial_temporal_hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    ode_dir = f"{dir}/{ode_model}/{config}/{closure_method}/closure_order_3/0.000000"
    save_dir = f"{save_dir}/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    plot_mean_errors(smm_dir, ode_dir, temporal_hybrid_dir, spatial_hybrid_dir, save_dir, tmax, comp_index, num_regions, figsize)
    plot_var_errors(smm_dir, ode_dir, temporal_hybrid_dir, spatial_hybrid_dir, save_dir, tmax, comp_index, num_regions, figsize)
