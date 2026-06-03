import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os
from settings import *
from cycler import cycler
import math

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
            ax.set_title(titles[r])
    fig.tight_layout()
    fig.savefig(f"{save_dir}/{list(compartment_colors.keys())[comp_index]}_sim_{start_sim}-{start_sim+num_runs}.png", dpi=dpi)
    plt.close(fig)

def plot_transitions(dir_path, save_dir, num_regions, figsize, tmax):
    os.makedirs(save_dir, exist_ok=True)

    # Preallocate arrays: (from, to, time)
    S = np.zeros((num_regions, num_regions, tmax))
    I = np.zeros((num_regions, num_regions, tmax))
    R = np.zeros((num_regions, num_regions, tmax))

    # --- Load data ---
    for t in range(tmax):
        matrix = pd.read_csv(f"{dir_path}/{t}_transitions.csv", header=None).values

        for fr in range(num_regions):
            for to in range(num_regions):
                if fr == to:
                    continue  # skip transition i->i

                S[fr, to, t] = matrix[fr * 3, to * 3]
                I[fr, to, t] = matrix[fr * 3 + 1, to * 3 + 1]
                R[fr, to, t] = matrix[fr * 3 + 2, to * 3 + 2]

    time = np.arange(tmax)

    # --- Helper to plot grid ---
    def plot_grid(data, filename, color):
        valid_pairs = [(fr, to) for fr in range(num_regions) for to in range(num_regions) if fr != to]
        n = len(valid_pairs)

        # Choose a roughly square layout
        ncols = math.ceil(np.sqrt(n))
        nrows = math.ceil(n / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.array(axes).reshape(-1)  # flatten safely

        for idx, (fr, to) in enumerate(valid_pairs):
            ax = axes[idx]
            ax.scatter(time, data[fr, to], color=color, s=10)
            ax.set_title(f"{fr}->{to}")

            # Left column → ylabel
            if idx % ncols == 0:
                ax.set_ylabel("Transitions [#]")

            # Bottom row → xlabel
            if idx >= (nrows - 1) * ncols:
                ax.set_xlabel("Time [days]")

        # Remove unused axes (if grid not perfect)
        for idx in range(len(valid_pairs), len(axes)):
            fig.delaxes(axes[idx])

        fig.tight_layout()
        fig.savefig(f"{save_dir}/{filename}")
        plt.close(fig)

    # --- Plot S, I, R grids ---
    plot_grid(S, "transitions_S.png", compartment_colors["S"][0])
    plot_grid(I, "transitions_I.png", compartment_colors["I"][0])
    plot_grid(R, "transitions_R.png", compartment_colors["R"][0])

    # --- Aggregated incoming transitions ---
    def plot_incoming(data, filename, color):
        # Sum over "from" axis → shape: (to, time)
        incoming = data.sum(axis=0)

        fig, axes = plt.subplots(num_regions, 1, figsize=(figsize[0] * 0.4, figsize[1] * 0.8))

        if num_regions == 1:
            axes = [axes]

        for to in range(num_regions):
            ax = axes[to]
            ax.plot(time, incoming[to], color=color)

            ax.set_title(f"Region {to}")

            ax.set_ylabel("Incoming [#]")

            if to == num_regions - 1:
                ax.set_xlabel("Time [days]")

        fig.tight_layout()
        fig.savefig(f"{save_dir}/{filename}")
        plt.close(fig)

    plot_incoming(S, "incoming_S.png", compartment_colors["S"][0])
    plot_incoming(I, "incoming_I.png", compartment_colors["I"][0])
    plot_incoming(R, "incoming_R.png", compartment_colors["R"][0])

def plot_mean(dir, save_dir, figsize, num_regions, comp_index):
    # Get number of figure rows and cols
    num_cols = 1
    i = 1
    while(i**2 < num_regions):
        num_cols += 1
        i += 1
    num_rows = min(num_cols, num_regions - num_cols + 1)
    # Create joint figure and axes
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    results = pd.read_csv(dir + f"/means.csv")
    for r in range(num_regions):
        if(num_regions == 1):
            ax = axes
        elif(num_rows == 1):
            ax = axes[r]
        else:
            ax = axes[int(r/num_rows), r%num_cols]
        
        ax.plot(results.Time, results.iloc[:, 1 + comp_index + r * len(compartment_names)], color = compartment_colors[list(compartment_colors.keys())[comp_index]][0])
        if(r%num_cols == 0):
            ax.set_ylabel(r"$\mu$")
        if(int(r/num_cols) == num_rows-1):
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])
        ax.set_title(titles[r])
    fig.tight_layout()
    fig.savefig(f"{save_dir}/{list(compartment_colors.keys())[comp_index]}_means.png", dpi=dpi)
    plt.close(fig)
    
def plot_std(dir, save_dir, figsize, num_regions, comp_index):
    # Get number of figure rows and cols
    num_cols = 1
    i = 1
    while(i**2 < num_regions):
        num_cols += 1
        i += 1
    num_rows = min(num_cols, num_regions - num_cols + 1)
    # Create joint figure and axes
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    results = pd.read_csv(dir + f"/moments.csv")
    for r in range(num_regions):
        if(num_regions == 1):
            ax = axes
        elif(num_rows == 1):
            ax = axes[r]
        else:
            ax = axes[int(r/num_rows), r%num_cols]
            
        col_name = "M"
        moment_index = [0 for _ in range(len(compartment_names) * num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        for m in moment_index:
            col_name += f"{m}"
        
        ax.plot(results.Time, np.sqrt(results[col_name].iloc[:]), color = compartment_colors[list(compartment_colors.keys())[comp_index]][0])
        if(r%num_cols == 0):
            ax.set_ylabel(r"$\sigma$")
        if(int(r/num_cols) == num_rows-1):
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])
        ax.set_title(titles[r])
    fig.tight_layout()
    fig.savefig(f"{save_dir}/{list(compartment_colors.keys())[comp_index]}_std.png", dpi=dpi)
    plt.close(fig)
            
if __name__ == "__main__":
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    model = "SMM"
    config = "config_SIR_I0_0_1_10_100_exchange"
    percentiles = ["05", "95"]
    other_model_dir = "V:/bick_ju/TemporalHybrid/Moments"
    other_label = "ODE"
    num_regions = 4
    tmax = 300
    figsize_percentiles = (4, 3)
    
    other_model_dir = f"{other_model_dir}/{config}/closure_order_2/0.000000/means.csv"
    dir = f"{dir}/{model}/{config}"
    save_dir = f"{save_dir}/{model}/{config}"
    os.makedirs(save_dir, exist_ok=True)
    titles = [r"Region 0, $I_0=0$, $R_0\approx4$", r"Region 1, $I_0=1$, $R_0\approx4$", r"Region 2, $I_0=10$, $R_0\approx4$", r"Region 3, $I_0=100$, $R_0\approx4$"]
    
    for i in range(0, 100, 10):
        #plot_simulation_results(10, dir, save_dir, (7, 5), num_regions, 0, i)
        plot_simulation_results(10, dir, save_dir, (7, 5), num_regions, 1, i)
    
    #plot_transitions(dir, save_dir, num_regions, (10, 9), tmax)
    #plot_mean(dir, save_dir, (7, 5), num_regions, 0)
    plot_mean(dir, save_dir, (7, 5), num_regions, 1)
    #plot_std(dir, save_dir, (7, 5), num_regions, 0)
    plot_std(dir, save_dir, (7, 5), num_regions, 1)

    # for r in range(num_regions):
    #     for comp_index in range(len(compartment_names)):
    #         plot_percentiles(dir, percentiles, save_dir, figsize_percentiles, comp_index + r * len(compartment_names), tmax, other_ts_file=other_model_dir, other_label=other_label)
