import pandas as pd
import matplotlib.pyplot as plt
import os
from settings import *

def plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index, figsize, tmax, num_regions):
    fig_var, ax_var = plt.subplots(figsize=figsize)
    fig_std, ax_std = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        moments = pd.read_csv(dir + f"/{config}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        time = moments['Time']
        ax_var.plot(time, moments[col_name], label=config_labels[c], color = compartment_colors[comp][c])
        ax_std.plot(time, np.sqrt(moments[col_name]), label=config_labels[c], color = compartment_colors[comp][c])

    ax_var.set_xlabel("Time [days]")
    ax_var.set_ylabel(f"Var {compartment_names[comp]}")
    ax_var.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax_var.legend()
    fig_var.tight_layout()
    fig_var.savefig(save_dir + f"/all_variances_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_var)
    
    ax_std.set_xlabel("Time [days]")
    ax_std.set_ylabel(f"Stddev {compartment_names[comp]}")
    ax_std.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax_std.legend()
    fig_std.tight_layout()
    fig_std.savefig(save_dir + f"/all_standard_deviations_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_std)

def plot_means(dir, configs, config_labels, save_dir, comp_index, figsize, tmax, num_regions):
    fig, ax = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    for c, config in enumerate(configs):        
        means = pd.read_csv(dir + f"/{config}/means.csv")
        means = means[means['Time'] <= tmax]
        time = means['Time']
        ax.plot(time, means.iloc[:, 1 + comp_index], label=config_labels[c], color = compartment_colors[comp][c])

    ax.set_xlabel("Time [days]")
    ax.set_ylabel(f"Mean {compartment_names[comp]}")
    ax.set_yscale('log')
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"/all_means_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig)

if __name__ == "__main__":
    figsize = (3.3, 2.5)#(4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    model = "Hybrid1"
    configs = ["config1_1r_I0_100/switch_value_0.000000", "config1_1r_I0_100/switch_value_0.010000", "config1_1r_I0_100/switch_value_1.000000"]
    config_labels = ["ODE", "Hybrid 1%", "SMM"]
    num_regions = 1
    tmax = 30
    
    dir = f"{dir}/{model}"
    save_dir = f"{save_dir}/{model}/{configs[0]}"
    os.makedirs(save_dir, exist_ok=True)
    
    plot_means(dir, configs, config_labels, save_dir, 1, figsize, tmax, num_regions)
    
    # for r in range(num_regions):
    #     for comp_index in range(len(compartment_names)):
    #         plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index + r * len(compartment_names), figsize, tmax, num_regions)
