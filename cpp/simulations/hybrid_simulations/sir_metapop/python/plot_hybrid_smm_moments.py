import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def plot_mean_ts(dir_hybrid, dir_smm, dir_mom, save_dir, comp_index, num_regions, figsize, tmax, tmin, color_hybrid, color_smm, color_mom):
    comp = list(compartment_colors.keys())[comp_index]
    
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        
        mean = pd.read_csv(dir_smm + f"/means.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_smm, label = r"Stochastic metapopulation model")
        
        # Plot time frame modeled with SMM
        mean = pd.read_csv(dir_hybrid + f"/{r}_means_smm.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_hybrid, label = r"Spatial-temporal-hybrid model")
        
        mean = pd.read_csv(dir_hybrid + f"/{r}_expected_values_moments.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index + r * len(compartment_names)], color=color_hybrid)
        # Add line at swithcing tp
        if(len(mean.Time) > 1):
            switch_tp = mean.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        mean = pd.read_csv(dir_mom + f"/means.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_mom, label = r"MoM")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\mu_I$")
        mean_smm = pd.read_csv(dir_smm + f"/means.csv")
        ax.set_ylim(-0.5*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)]), 2*np.max(mean_smm.iloc[1:, 1 + comp_index]))
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        
        handles, labels = ax.get_legend_handles_labels()
        fig_leg = plt.figure(figsize=figsize)
        fig_leg.legend(handles, labels, loc='center')
        fig_leg.tight_layout()
        fig_leg.savefig(f"{save_dir}/legend_mean_{compartment_names[comp]}_r{r}.png", dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)


def plot_var_ts(dir_hybrid, dir_smm, dir_mom, save_dir, comp_index, num_regions, figsize, tmax, tmin, color_hybrid, color_smm, color_mom):
    comp = list(compartment_colors.keys())[comp_index]
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
            
        fig, ax = plt.subplots(figsize=figsize)        
        var = pd.read_csv(dir_smm + f"/moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color=color_smm, label = r"SMM")
        
        var = pd.read_csv(dir_hybrid + f"/{r}_moments_smm.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color = color_hybrid, label = r"STHMM")
        
        var = pd.read_csv(dir_hybrid + f"/{r}_moments_moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = color_hybrid)
        # Add line at swithcing tp
        if(len(var.Time) > 1):
            switch_tp = var.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        var = pd.read_csv(dir_mom + f"/moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color=color_mom, label = r"MoM")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\sigma^2_I$")
        var_smm = pd.read_csv(dir_smm + f"/moments.csv")
        ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
        #ax.set_yscale("log")
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)

def plot_runtimes(dir_hybrid, dir_smm, dir_mom, save_dir, figsize, color_hybrid, color_smm, color_mom):
    times = []
    time = pd.read_csv(dir_smm + f"/total_time.csv")
    times.append(time["Runtime"].iloc[0])
    time = pd.read_csv(dir_hybrid + f"/total_time.csv")
    times.append(time["Runtime"].iloc[0])
    time = pd.read_csv(dir_mom + f"/total_time.csv")
    times.append(time["Runtime"].iloc[0])
    
    print("Speed up hybrid: " + str(times[0]/times[1]), "Gain: " + str((times[0] - times[1])/times[0] * 100) + "%")
    
    fig, ax = plt.subplots(figsize=figsize)
    x = [0, 1, 2]
    bar_colors = [color_smm, color_hybrid, color_mom]
    ax.bar(x, times, color =  bar_colors, width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Runtime [s]")
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.1, top=0.95)
    fig.savefig(f"{save_dir}/runtimes.png", dpi=dpi)
    plt.close(fig)
    

if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Spatial-Hybrid2"
    config = "config_2r_100_0_k1"
    num_regions = 2
    condition = "var_gradient"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    color_hybrid = colors['middle green']
    color_mom = colors['brown']
    tmin = 0
    tmax = 200
    condition_name = ""
    
    smm_dir = f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    mom_dir = f"{dir}/Moments/{config}/{closure_method}/{closure_order}/0.000000"
    save_dir = f"{save_dir}/Comparison/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    # plot_mean_ts(hybrid_dir, smm_dir, mom_dir, save_dir, comp_index=1, num_regions=num_regions, figsize=figsize, tmax=tmax, tmin=tmin, color_hybrid=color_hybrid, color_smm=color_smm, color_mom=color_mom)
    # plot_var_ts(hybrid_dir, smm_dir, mom_dir, save_dir, comp_index=1, num_regions=num_regions, figsize=figsize, tmax=tmax, tmin=tmin, color_hybrid=color_hybrid, color_smm=color_smm, color_mom=color_mom)
    
    plot_runtimes(hybrid_dir, smm_dir, mom_dir, save_dir, figsize=figsize, color_hybrid=color_hybrid, color_smm=color_smm, color_mom=color_mom)
