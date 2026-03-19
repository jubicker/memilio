import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def plot_mean_ts(dir, save_dir, comp_index, num_regions, figsize, tmax, tmin, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_smm, label = r"SMM")
        
        # Plot time frame modeled with SMM
        mean = pd.read_csv(dir + f"/{r}_means_smm.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=compartment_colors[comp][0], label = r"STHMM")
        
        mean = pd.read_csv(dir + f"/{r}_expected_values_moments.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index + r * len(compartment_names)], color=compartment_colors[comp][0])
        # Add line at swithcing tp
        if(len(mean.Time) > 1):
            switch_tp = mean.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\mu_I$")
        if(dir_smm!=""):
            mean_smm = pd.read_csv(dir_smm + f"/means.csv")
            ax.set_ylim(-0.5*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)]), 2*np.max(mean_smm.iloc[1:, 1 + comp_index]))
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
    
    if(num_regions > 1):
        # One plot for both regions
        fig, ax = plt.subplots(figsize=figsize)
        colors_smm = [colors["very dark blue"], colors["very dark green"]]
        colors_hybrid = [colors["middle blue"], colors["middle green"]]
        for r in range(num_regions):
                      
            # Plot time frame modeled with SMM
            mean = pd.read_csv(dir + f"/{r}_means_smm.csv")
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=colors_hybrid[r], label = f"STHMM region {r}")
            
            mean = pd.read_csv(dir + f"/{r}_expected_values_moments.csv")
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index + r * len(compartment_names)], color=colors_hybrid[r])
            
            if(dir_smm) != "":
                mean = pd.read_csv(dir_smm + f"/means.csv")
                mean = mean[mean['Time'] <= tmax]
                mean = mean[mean['Time'] >= tmin]
                ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=colors_smm[r], label = f"SMM region {r}", linestyle="dotted")
                
            ax.set_xlabel("Time [days]")
            ax.set_ylabel(r"$\mu_I$")
            if(dir_smm!=""):
                mean_smm = pd.read_csv(dir_smm + f"/means.csv")
                ax.set_ylim(-0.1*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)]), 1.1*np.max(mean_smm.iloc[1:, 1 + comp_index]))
            ax.legend()
            fig.tight_layout()
            fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_all_regions.png", dpi=dpi)
            plt.close(fig)
            
    # One plot for both regions
    fig, ax = plt.subplots(figsize=figsize)
    region_colors_hybrid = [colors["middle blue"], colors["middle green"]]
    region_colors_smm = [colors["dark blue"], colors["dark green"]]
    for r in range(num_regions):
        
        # Plot time frame modeled with SMM
        mean = pd.read_csv(dir + f"/{r}_means_smm.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=region_colors_hybrid[r], label = r"Region " + str(r+1), linewidth=2)
        
        mean = pd.read_csv(dir + f"/{r}_expected_values_moments.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index + r * len(compartment_names)], color=region_colors_hybrid[r])
        
        # Add line at swithcing tp
        if(len(mean.Time) > 1):
            switch_tp = mean.Time.iloc[1]
            ax.axvline(x=switch_tp, color=region_colors_smm[r], linestyle="--")
            
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=region_colors_smm[r], label = r"Stochastic model, region " + str(r+1), linestyle="dotted", linewidth=2)
            
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(r"$\mu_I$")
    #ax.legend()
    fig.tight_layout()
    fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_all_regions.png", dpi=dpi)
    plt.close(fig)
        
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/all_regions_legend.png", dpi=dpi, bbox_inches='tight')
        
def plot_var_ts(dir, save_dir, comp_index, num_regions, figsize, tmax, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
            
        fig, ax = plt.subplots(figsize=figsize)        
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            var = var[var['Time'] <= tmax]
            ax.plot(var.Time, var[col_name].iloc[:], color=color_smm, label = r"SMM")
        
        # Plot SMM variance
        var = pd.read_csv(dir + f"/{r}_moments_smm.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color = compartment_colors[comp][0], label = r"STHMM")
        
        # Plot moment variance
        var = pd.read_csv(dir + f"/{r}_moments_moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = compartment_colors[comp][0])
        # Add line at swithcing tp
        if(len(var.Time) > 1):
            switch_tp = var.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\sigma^2_I$")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
        #ax.set_yscale("log")
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
        
    # All regions in one plot
    region_colors_hybrid = [colors["middle blue"], colors["middle green"]]
    region_colors_smm = [colors["dark blue"], colors["dark green"]]
    fig, ax = plt.subplots(figsize=figsize)    
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        
        # Plot SMM variance
        var = pd.read_csv(dir + f"/{r}_moments_smm.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color = region_colors_hybrid[r], label = r"THMM, region " + str(r+1), linewidth=2)
        
        # Plot moment variance
        var = pd.read_csv(dir + f"/{r}_moments_moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = region_colors_hybrid[r], linewidth=2)
        # Add line at swithcing tp
        if(len(var.Time) > 1):
            switch_tp = var.Time.iloc[1]
            ax.axvline(x=switch_tp, color=region_colors_smm[r], linestyle="--")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\sigma^2_I$")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
            
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            var = var[var['Time'] <= tmax]
            ax.plot(var.Time, var[col_name].iloc[:], color=region_colors_smm[r], label = r"SMM, region " + str(r+1), linestyle="dotted", linewidth=2)
        #ax.set_yscale("log")
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_all_regions.png", dpi=dpi)
        plt.close(fig)


if __name__ == "__main__":
    figsize = (4.5, 3.2)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Spatial-Hybrid2"
    config = "config_2r_10_100_k1"
    num_regions = 2
    condition = "var_gradient"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    tmin = 0
    tmax = 90
    condition_name = ""
    
    smm_dir = ""#f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    plot_mean_ts(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax, tmin, smm_dir, color_smm)
    plot_var_ts(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax, smm_dir, color_smm)
