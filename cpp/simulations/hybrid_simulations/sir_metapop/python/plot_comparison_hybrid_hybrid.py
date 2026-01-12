import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def plot_mean_error(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir, save_dir, switching_thresholds, colors_hybrid1, colors_hybrid2, color_ode, comp_index, figsize):
    comp = list(compartment_colors.keys())[comp_index]
    fig, ax = plt.subplots(figsize=figsize)
    
    smm_values = pd.read_csv(smm_dir + "/means.csv")
    ode_values = pd.read_csv(ode_dir + "/means.csv")
    
    summed_errors = []
    
    error = (smm_values.iloc[:-1, 1 + comp_index]-ode_values.iloc[:, 1 + comp_index]).abs()
    summed_errors.append(error.sum())
    
    ax.plot(ode_values.Time, error, label = "ODE", color = color_ode)
    
    for i, st in enumerate(switching_thresholds):
        values = pd.read_csv(hybrid1_dir + f"/switch_value_{st}/means.csv")
        error = (smm_values.iloc[:-1, 1 + comp_index] - values.iloc[:, 1 + comp_index]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.2f}% Hybrid ODE", color = colors_hybrid1[i])
        
    for i, st in enumerate(switching_thresholds):
        values = pd.read_csv(hybrid2_dir + f"/switch_value_{st}/means.csv")
        error = (smm_values.iloc[:-1, 1 + comp_index] - values.iloc[:, 1 + comp_index]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.2f} Hybrid MoM%", color = colors_hybrid2[i])
        
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
    x = range(2*len(switching_thresholds) + 1)
    bar_colors = [color_ode] + colors_hybrid1 + colors_hybrid2
    ax.bar(x, summed_errors, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/mean_sum_err_{comp}.png", dpi=dpi)
    plt.close()

def plot_variance_error(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir ,save_dir, switching_thresholds, colors_hybrid1, colors_hybrid2, color_ode, comp_index, figsize):
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
        values = pd.read_csv(hybrid1_dir + f"/switch_value_{st}/moments.csv")
        error = (smm_values[col_name].iloc[:-1] - values[col_name].iloc[:]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.2f}% Hybrid ODE", color = colors_hybrid1[i])
        
    for i, st in enumerate(switching_thresholds):
        values = pd.read_csv(hybrid2_dir + f"/switch_value_{st}/moments.csv")
        error = (smm_values[col_name].iloc[:-1] - values[col_name].iloc[:]).abs()
        summed_errors.append(error.sum())
        ax.plot(values.Time, error, label = f"{float(st) * 100:.2f}% Hybrid MoM", color = colors_hybrid2[i])
        
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
    x = range(2*len(switching_thresholds) + 1)
    bar_colors = [color_ode] + colors_hybrid1 + colors_hybrid2
    ax.bar(x, summed_errors, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/var_sum_err_{comp}.png", dpi=dpi)
    plt.close()

def plot_runtimes(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir, save_dir, switching_thresholds, colors_hybrid1, colors_hybrid2, color_ode, color_smm, figsize):    
    total_runtime = []
    labels = []
    
    # ODE
    runtime = pd.read_csv(ode_dir + "/total_time.csv")
    total_runtime.append(runtime['Runtime'].iloc[0])
    labels.append("ODE")
    # Hybrid1
    for st in switching_thresholds:
        runtime = pd.read_csv(hybrid1_dir + f"/switch_value_{st}/total_time.csv")
        total_runtime.append(runtime['Runtime'].iloc[0])
        labels.append(f"{float(st) * 100:.2f}% Hybrid ODE")
        
    for st in switching_thresholds:
        runtime = pd.read_csv(hybrid2_dir + f"/switch_value_{st}/total_time.csv")
        total_runtime.append(runtime['Runtime'].iloc[0])
        labels.append(f"{float(st) * 100:.2f}% Hybrid MoM")
        
    # SMM    
    runtime = pd.read_csv(smm_dir + "/total_time.csv")
    total_runtime.append(runtime['Runtime'].iloc[0])
    labels.append("SMM")
    
    bar_colors = [color_ode] + colors_hybrid1 + colors_hybrid2 + [color_smm]
    
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(2*len(switching_thresholds) + 2)
    ax.bar(x, total_runtime, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Runtime [s]")
    #ax.set_yscale("log")
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


if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel/Hybrid_comparison"
    hybrid1_model = "Hybrid1"
    hybrid2_model = "Hybrid2"
    config = "config1_1r_I0_1"
    switching_values = ["0.001000", "0.010000", "0.100000", "0.300000"]
    colors_hybrid1 = [colors['dark blue'], colors['middle blue'], colors['light blue'], colors['teal']]
    colors_hybrid2 = [colors['purple'], colors['rose'], colors['red'], colors['dark red']]
    color_ode = "black"
    color_smm = colors['dark grey']
    
    smm_dir = f"{dir}/SMM/{config}"
    hybrid1_dir = f"{dir}/{hybrid1_model}/{config}"
    hybrid2_dir = f"{dir}/{hybrid2_model}/{config}"
    ode_dir = f"{dir}/Moments/{config}/closure_order_2/0.000000"
    save_dir = f"{save_dir}/{config}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    
    plot_mean_error(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir, save_dir, switching_values, colors_hybrid1, colors_hybrid2, color_ode, comp_index, figsize)
    plot_variance_error(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir, save_dir, switching_values, colors_hybrid1, colors_hybrid2, color_ode, comp_index, figsize)
    plot_runtimes(smm_dir, hybrid1_dir, hybrid2_dir, ode_dir, save_dir, switching_values, colors_hybrid1, colors_hybrid2, color_ode, color_smm, figsize)
