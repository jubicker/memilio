import pandas as pd
import matplotlib.pyplot as plt
import os
from settings import *

def plot_percentiles(dir, percentiles, save_dir, figsize, comp_index, other_ts_file = "", other_label="", color_other = "black"):
       
    fig, ax = plt.subplots(figsize=figsize)
    alpha = 0.3
    comp = list(compartment_colors.keys())[comp_index]
    
    # Plot percentiles
    while(len(percentiles) > 0):
        lower = pd.read_csv(f"{dir}/p{percentiles[0]}.csv")
        upper = pd.read_csv(f"{dir}/p{percentiles[-1]}.csv")
        time = lower['Time']
        label = int(percentiles[-1]) - int(percentiles[0])
        ax.fill_between(time, lower.iloc[:, comp_index + 1], upper.iloc[:, comp_index + 1], color=compartment_colors[comp][0], alpha=alpha, label=f"{label}% percentile")
        
        # Remove used percentiles
        percentiles = percentiles[1:-1]
        
        alpha += 0.2
        
    # Plot mean
    mean = pd.read_csv(f"{dir}/means.csv")
    ax.plot(mean['Time'], mean.iloc[:, comp_index + 1], color=compartment_colors[comp][0], label="Mean", linewidth=1)
    
    # Plot other ts if provided
    if other_ts_file != "":
        other_ts = pd.read_csv(other_ts_file)
        ax.plot(other_ts['Time'], other_ts.iloc[:, comp_index + 1], color=color_other, linestyle="--", label=other_label, linewidth=1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(f"{compartment_names[comp]} [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{save_dir}/percentiles_{compartment_names[comp]}.png", dpi=dpi)
    plt.close(fig)
    
if __name__ == "__main__":
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    model = "SMM"
    config = "config1_1r_I0_2"
    percentiles = ["05", "95"]
    other_model_dir = "V:/bick_ju/TemporalHybrid/Moments"
    other_label = "ODE"
    figsize_percentiles = (3.3, 2.5)
    
    other_model_dir = f"{other_model_dir}/{config}/closure_order_2/0.000000/means.csv"
    dir = f"{dir}/{model}/{config}"
    save_dir = f"{save_dir}/{model}/{config}"
    os.makedirs(save_dir, exist_ok=True)

    for comp_index in range(len(compartment_names)):
        plot_percentiles(dir, percentiles, save_dir, figsize_percentiles, comp_index, other_ts_file=other_model_dir, other_label=other_label)
