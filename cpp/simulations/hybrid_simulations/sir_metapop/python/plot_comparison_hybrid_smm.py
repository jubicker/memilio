import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def plot_mean_error(smm_dir, hybrid_dir, ode_dir,save_dir, switching_thresholds, colors_hybrid, color_ode, comp_index, figsize):
    comp = list(compartment_colors.keys())[comp_index]
    fig, ax = plt.subplots(figsize=figsize)
    
    smm_values = pd.read_csv(smm_dir + "/means.csv")
    ode_values = pd.read_csv(ode_dir + "/means.csv")
    
    summed_errors = []
    
    error = (smm_values.iloc[:-1, 1 + comp_index]-ode_values.iloc[:, 1 + comp_index]).abs()
    summed_errors.append(error.sum())
    
    ax.plot(ode_values.Time, error, label = "ODE", color = color_ode)
    
    for i, st in enumerate(switching_thresholds):
        values = pd.read_csv(hybrid_dir + f"/switch_value_{st}/means.csv")
        error = (smm_values.iloc[:-1, 1 + comp_index] - values.iloc[:, 1 + comp_index]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.4f}%", color = colors_hybrid[i])
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Log(err)")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/mean_error_ts_{comp}")
    
    # Save legend
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center')      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + f"/legend_{comp}.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    # Save summed error
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(switching_thresholds) + 1)
    bar_colors = [color_ode] + colors_hybrid
    ax.bar(x, summed_errors, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/mean_sum_err_{comp}.png", dpi=dpi)
    plt.close()
    
def plot_variance_error(smm_dir, hybrid_dir, ode_dir,save_dir, switching_thresholds, colors_hybrid, color_ode, comp_index, figsize):
    comp = list(compartment_colors.keys())[comp_index]
    moment_index = [0 for _ in range(len(compartment_names))]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M{moment_index[0]}{moment_index[1]}{moment_index[2]}"
    fig, ax = plt.subplots(figsize=figsize)
    
    smm_values = pd.read_csv(smm_dir + "/moments.csv")
    ode_values = pd.read_csv(ode_dir + "/moments.csv")
    
    summed_errors = []
    
    error = (smm_values[col_name].iloc[:-1] - ode_values[col_name].iloc[:]).abs()
    summed_errors.append(error.sum())
    
    ax.plot(ode_values.Time, error, label = "ODE", color = color_ode)
    
    for i, st in enumerate(switching_thresholds):
        values = pd.read_csv(hybrid_dir + f"/switch_value_{st}/moments.csv")
        error = (smm_values[col_name].iloc[:-1] - values[col_name].iloc[:]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.4f}%", color = colors_hybrid[i])
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Log(err)")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/var_error_ts_{comp}")
    
    # Save legend
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center')      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + f"/legend_{comp}.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    # Save summed error
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(switching_thresholds) + 1)
    bar_colors = [color_ode] + colors_hybrid
    ax.bar(x, summed_errors, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/var_sum_err_{comp}.png", dpi=dpi)
    plt.close()

def plot_runtimes(smm_dir, hybrid_dir, ode_dir, save_dir, switching_thresholds, colors_hybrid, color_ode, color_smm, figsize):    
    mean_runtimes = []
    max_runtimes = []
    total_runtime = []
    labels = []
    
    # ODE
    runtime = pd.read_csv(ode_dir + "/total_time.csv")
    mean_runtimes.append(runtime['Runtime'].iloc[0])
    max_runtimes.append(runtime['Runtime'].iloc[0])
    total_runtime.append(runtime['Runtime'].iloc[0])
    labels.append("ODE")
    # Hybrid
    for st in switching_thresholds:
        runtime = pd.read_csv(hybrid_dir + f"/switch_value_{st}/runtimes.csv")
        mean_runtimes.append(runtime['Runtime'].mean())
        max_runtimes.append(runtime['Runtime'].max())
        runtime = pd.read_csv(hybrid_dir + f"/switch_value_{st}/total_time.csv")
        total_runtime.append(runtime['Runtime'].iloc[0])
        labels.append(f"{float(st) * 100:.4f}%")
        
    # SMM    
    runtime  = pd.read_csv(smm_dir + "/runtimes.csv")
    mean_runtimes.append(runtime['Runtime'].mean())
    max_runtimes.append(runtime['Runtime'].max())
    runtime = pd.read_csv(smm_dir + "/total_time.csv")
    total_runtime.append(runtime['Runtime'].iloc[0])
    labels.append("SMM")
    
    bar_colors = [color_ode] + colors_hybrid + [color_smm]
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(switching_thresholds) + 2)
    ax.bar(x, mean_runtimes, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/mean_runtime.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(switching_thresholds) + 2)
    ax.bar(x, max_runtimes, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Runtime [s]")
    #ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/max_runtime.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(switching_thresholds) + 2)
    ax.bar(x, total_runtime, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/total_runtime.png", dpi=dpi)
    plt.close()
    
    # Legend
    handles = [Patch(color=c, label=l) for c, l in zip(bar_colors, labels)]
    fig, ax = plt.subplots(figsize=figsize)
    ax.legend(handles=handles, loc='center')
    ax.axis('off')
    fig.savefig(save_dir + '/legend_runtimes.png', bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    
def plot_mean_var_ts(dir, name_mean, name_var, switching_thresholds, comp_index, figsize, name2_mean="", name2_var="", dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    moment_index = [0 for _ in range(len(compartment_names))]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M{moment_index[0]}{moment_index[1]}{moment_index[2]}"
    for st in switching_thresholds:
        fig, ax = plt.subplots(figsize=figsize)
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            var = pd.read_csv(dir_smm + f"/moments.csv")
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index], color=color_smm, label = r"$\mu$ SMM")
            ax.fill_between(mean.Time, mean.iloc[:, 1 + comp_index] - np.sqrt(var[col_name].iloc[:]), mean.iloc[:, 1 + comp_index] + np.sqrt(var[col_name].iloc[:]), color = color_smm, alpha = 0.5)
            ax.fill_between(mean.Time, mean.iloc[:, 1 + comp_index] - 2*np.sqrt(var[col_name].iloc[:]), mean.iloc[:, 1 + comp_index] + 2*np.sqrt(var[col_name].iloc[:]), color = color_smm, alpha = 0.2)
            
        mean = pd.read_csv(dir + f"/switch_value_{st}/" + name_mean + ".csv")
        var = pd.read_csv(dir + f"/switch_value_{st}/" + name_var + ".csv")
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index], color=compartment_colors[comp][0], label = r"$\mu$")
        ax.fill_between(mean.Time, mean.iloc[:, 1 + comp_index] - np.sqrt(var[col_name].iloc[:]), mean.iloc[:, 1 + comp_index] + np.sqrt(var[col_name].iloc[:]), color = compartment_colors[comp][0], alpha = 0.5, label = r"$\mu\pm\sigma$")
        ax.fill_between(mean.Time, mean.iloc[:, 1 + comp_index] - 2*np.sqrt(var[col_name].iloc[:]), mean.iloc[:, 1 + comp_index] + 2*np.sqrt(var[col_name].iloc[:]), color = compartment_colors[comp][0], alpha = 0.2, label = r"$\mu\pm 2\sigma$")
        
        if(name2_mean != ""):
            mean = pd.read_csv(dir + f"/switch_value_{st}/" + name2_mean + ".csv")
            var = pd.read_csv(dir + f"/switch_value_{st}/" + name2_var + ".csv")
            ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index], color=compartment_colors[comp][0])
            ax.fill_between(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index] - np.sqrt(var[col_name].iloc[1:]), mean.iloc[1:, 1 + comp_index] + np.sqrt(var[col_name].iloc[1:]), color = compartment_colors[comp][0], alpha = 0.5)
            ax.fill_between(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index] - 2*np.sqrt(var[col_name].iloc[1:]), mean.iloc[1:, 1 + comp_index] + 2*np.sqrt(var[col_name].iloc[1:]), color = compartment_colors[comp][0], alpha = 0.2)
            # Add line at swithcing tp
            switch_tp = mean.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
        
        if(dir_smm!=""):
            mean_smm = pd.read_csv(dir_smm + f"/means.csv")    
            ax.set_ylim(-0.01*np.max(mean_smm.iloc[:, 1 + comp_index]), 1.5*np.max(mean_smm.iloc[:, 1 + comp_index]))
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(f"{compartment_names[comp]} [#]")
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_var_{compartment_names[comp]}_{st}.png", dpi=dpi)
        plt.close(fig)
        
def plot_mean_ts(dir, name_mean, switching_thresholds, comp_index, figsize, name2_mean="", dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    moment_index = [0 for _ in range(len(compartment_names))]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    for st in switching_thresholds:
        fig, ax = plt.subplots(figsize=figsize)
        
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index], color=color_smm, label = r"$\mu$ SMM")
        
        mean = pd.read_csv(dir + f"/switch_value_{st}/" + name_mean + ".csv")
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index], color=compartment_colors[comp][0], label = r"$\mu$")
        
        if(name2_mean != ""):
            mean = pd.read_csv(dir + f"/switch_value_{st}/" + name2_mean + ".csv")
            ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index], color=compartment_colors[comp][0])
            # Add line at swithcing tp
            switch_tp = mean.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(f"{compartment_names[comp]} [#]")
        if(dir_smm!=""):
            mean_smm = pd.read_csv(dir_smm + f"/means.csv")
            ax.set_ylim(-0.5*np.max(mean_smm.iloc[1:, 1 + comp_index]), 2*np.max(mean_smm.iloc[1:, 1 + comp_index]))
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_{st}.png", dpi=dpi)
        plt.close(fig)
        
def plot_var_ts(dir, name_var, switching_thresholds, comp_index, figsize, name2_var="", dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    moment_index = [0 for _ in range(len(compartment_names))]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M{moment_index[0]}{moment_index[1]}{moment_index[2]}"
    for st in switching_thresholds:
        fig, ax = plt.subplots(figsize=figsize)
        
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            ax.plot(var.Time, var[col_name].iloc[:], color=color_smm, label = r"$\sigma^2$ SMM")
        
        var = pd.read_csv(dir + f"/switch_value_{st}/" + name_var + ".csv")
        ax.plot(var.Time, var[col_name].iloc[:], color = compartment_colors[comp][0], label = r"$\sigma^2$")
        
        if(name2_var != ""):
            var = pd.read_csv(dir + f"/switch_value_{st}/" + name2_var + ".csv")
            ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = compartment_colors[comp][0])
            # Add line at swithcing tp
            switch_tp = var.Time.iloc[1]
            ax.axvline(x=switch_tp, color="black", linestyle="--")
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(f"{compartment_names[comp]} [#]")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
        #ax.set_yscale("log")
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_{st}.png", dpi=dpi)
        plt.close(fig)
     

if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Hybrid2"
    config = "config1_1r_I0_1"
    switching_values = ["0.000000", "0.000010", "0.000100", "0.001000", "0.010000", "0.100000", "0.300000"]
    colors_hybrid = [colors['dark blue'], colors['light blue'], colors['teal'], colors['light teal'], colors['dark green'], colors['middle green'], colors['light green'], colors['dark grey'], colors['middle grey']]
    color_ode = "black"
    color_smm = colors['dark grey']
    
    smm_dir = f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}"
    ode_dir = f"{dir}/Moments/{config}/closure_order_2/0.000000"
    save_dir = f"{save_dir}/{hybrid_model}/{config}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    
    plot_mean_ts(hybrid_dir, "means_smm", ["1.000000"], comp_index, figsize)
    plot_var_ts(hybrid_dir, "moments_smm", ["1.000000"], comp_index, figsize)
    plot_mean_var_ts(hybrid_dir, "means_smm", "moments_smm", ["1.000000"], comp_index, figsize)
    plot_mean_ts(hybrid_dir, "means_smm", switching_values, comp_index, figsize, "expected_values_moments", smm_dir, color_smm)
    plot_var_ts(hybrid_dir, "moments_smm", switching_values, comp_index, figsize, "moments_moments", smm_dir, color_smm)
    plot_mean_var_ts(hybrid_dir, "means_smm", "moments_smm", switching_values, comp_index, figsize, "expected_values_moments", "moments_moments", smm_dir, color_smm)
    
    #for comp_index in range(len(compartment_names)):
    plot_mean_error(smm_dir, hybrid_dir, ode_dir, save_dir, switching_values, colors_hybrid, color_ode, comp_index, figsize)
    plot_variance_error(smm_dir, hybrid_dir, ode_dir, save_dir, switching_values, colors_hybrid, color_ode, comp_index, figsize)
    plot_runtimes(smm_dir, hybrid_dir, ode_dir, save_dir, switching_values, colors_hybrid, color_ode, color_smm, figsize)
