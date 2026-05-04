import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os
from settings import *
from cycler import cycler

def plot_percentiles(dir, percentiles, save_dir, figsize, comp_index, tmax, other_ts_file = "", other_label="", color_other = "black"):
       
    fig, ax = plt.subplots(figsize=figsize)
    alpha = 0.3
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    
    # Plot percentiles
    while(len(percentiles) > 0):
        lower = pd.read_csv(f"{dir}/p{percentiles[0]}.csv")
        lower = lower[lower['Time'] <= tmax]
        upper = pd.read_csv(f"{dir}/p{percentiles[-1]}.csv")
        upper = upper[upper['Time'] <= tmax]
        time = lower['Time']
        label = int(percentiles[-1]) - int(percentiles[0])
        ax.fill_between(time, lower.iloc[:, comp_index + 1], upper.iloc[:, comp_index + 1], color=compartment_colors[comp][0], alpha=alpha, label=f"{label}% percentile")
        
        # Remove used percentiles
        percentiles = percentiles[1:-1]
        
        alpha += 0.2
        
    # Plot mean
    mean = pd.read_csv(f"{dir}/means.csv")
    mean = mean[mean['Time'] <= tmax]
    ax.plot(mean['Time'], mean.iloc[:, comp_index + 1], color=compartment_colors[comp][0], label="Mean", linewidth=1)
    
    # Plot other ts if provided
    if other_ts_file != "":
        other_ts = pd.read_csv(other_ts_file)
        other_ts = other_ts[other_ts['Time'] <= tmax]
        ax.plot(other_ts['Time'], other_ts.iloc[:, comp_index + 1], color=color_other, linestyle="--", label=other_label, linewidth=1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(f"{compartment_names[comp]} [#]")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    fig.tight_layout()
    fig.savefig(f"{save_dir}/percentiles_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    
    _, labels = ax.get_legend_handles_labels()
    proxy_handles = []
    for label in labels:
        if "percentile" in label.lower():
            proxy_handles.append(
                Patch(
                    facecolor=compartment_colors[comp][0],
                    alpha=0.5,
                    edgecolor="none"
                )
            )
        elif label == other_label:
            proxy_handles.append(
                Line2D(
                    [0], [0],
                    color=color_other,
                    linestyle="--",
                    linewidth=1
                )
            )
        else:
            proxy_handles.append(
                Line2D(
                    [0], [0],
                    color=compartment_colors[comp][0],
                    linewidth=1
                )
            )
    fig_leg, ax_leg = plt.subplots(figsize=(4, 3))
    ax_leg.axis("off")
    ax_leg.legend(proxy_handles, labels, loc="center")
    fig_leg.savefig(
        f"{save_dir}/legend_{compartment_names[comp]}_r{region}.png",
        dpi=dpi,
        bbox_inches="tight"
    )
    plt.close(fig)
    plt.close(fig_leg)
 
def plot_simulation_results(num_runs, dir, save_dir, figsize, num_regions, comp_index, start_sim=0):
    # Get number of figure rows and cols
    num_cols = 1
    i = 1
    while(i**2 < num_regions):
        num_cols += 1
        i += 1
    num_rows = min(num_cols, num_regions - num_cols + 1)
    # Create joint figure and axes
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    for sim in range(start_sim, start_sim + num_runs):
        results = pd.read_csv(dir + f"/{sim}_result.csv")
        for r in range(num_regions):
            if(num_regions == 1):
                ax = axes
            elif(num_rows == 1):
                ax = axes[r]
            else:
                ax = axes[int(r/num_rows), r%num_cols]
            
            ax.plot(results.Time, results.iloc[:, 1 + comp_index + r * len(compartment_names)])
            if(r%num_cols == 0):
                ax.set_ylabel(f"{compartment_names[list(compartment_colors.keys())[comp_index]]} [#]")
            if(int(r/num_cols) == num_rows-1):
                ax.set_xlabel("Time [days]")
            else:
                ax.set_xticks([])
            ax.set_title(f"Region {r}")
    fig.tight_layout()
    fig.savefig(f"{save_dir}/{list(compartment_colors.keys())[comp_index]}_sim_{start_sim}-{start_sim+num_runs}.png", dpi=dpi)
    plt.close(fig)
        
if __name__ == "__main__":
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    model = "SMM"
    config = "SIR/config_4r_0_1_10_100_k2"
    percentiles = ["05", "95"]
    other_model_dir = "V:/bick_ju/TemporalHybrid/Moments"
    other_label = "ODE"
    num_regions = 2
    tmax = 250
    figsize_percentiles = (4, 3)
    
    other_model_dir = f"{other_model_dir}/{config}/closure_order_2/0.000000/means.csv"
    dir = f"{dir}/{model}/{config}"
    save_dir = f"{save_dir}/{model}/{config}"
    os.makedirs(save_dir, exist_ok=True)
    
    plot_simulation_results(100, dir, save_dir, (7, 5), 4, 1, 0)

    # for r in range(num_regions):
    #     for comp_index in range(len(compartment_names)):
    #         plot_percentiles(dir, percentiles, save_dir, figsize_percentiles, comp_index + r * len(compartment_names), tmax, other_ts_file=other_model_dir, other_label=other_label)
