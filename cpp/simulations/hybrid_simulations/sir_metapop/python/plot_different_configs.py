import pandas as pd
import matplotlib.pyplot as plt
import os
from settings import *

def plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index, figsize):
    fig_var, ax_var = plt.subplots(figsize=figsize)
    fig_std, ax_std = plt.subplots(figsize=figsize)
    comp = list(compartment_names.keys())[comp_index]
    moment_index = [0 for _ in range(len(compartment_names))]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M{moment_index[0]}{moment_index[1]}{moment_index[2]}"
    for c, config in enumerate(configs):        
        moments = pd.read_csv(dir + f"/{config}/moments.csv")
        time = moments['Time']
        ax_var.plot(time, moments[col_name], label=config_labels[c], color = compartment_colors[comp][c])
        ax_std.plot(time, np.sqrt(moments[col_name]), label=config_labels[c], color = compartment_colors[comp][c])

    ax_var.set_xlabel("Time [days]")
    ax_var.set_ylabel(f"Var {compartment_names[comp]}")
    ax_var.legend()
    fig_var.tight_layout()
    fig_var.savefig(save_dir + f"/all_variances_{compartment_names[comp]}.png", dpi=dpi)
    plt.close(fig_var)
    
    ax_std.set_xlabel("Time [days]")
    ax_std.set_ylabel(f"Stddev {compartment_names[comp]}")
    ax_std.legend()
    fig_std.tight_layout()
    fig_std.savefig(save_dir + f"/all_standard_deviations_{compartment_names[comp]}.png", dpi=dpi)
    plt.close(fig_std)

if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    model = "SMM"
    configs = ["config1_1r_I0_1", "config1_1r_I0_2", "config1_1r_I0_10", "config1_1r_I0_100"]
    config_labels = ["I0=1", "I0=2", "I0=10", "I0=100"]
    
    dir = f"{dir}/{model}"
    save_dir = f"{save_dir}/{model}"
    os.makedirs(save_dir, exist_ok=True)
    
    for comp_index in range(len(compartment_names)):
        plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index, figsize)
