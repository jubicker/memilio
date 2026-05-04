import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def get_shaded_intervals(mask_time, mask_vals):
    intervals = []
    start = None

    for i in range(len(mask_vals)):
        if mask_vals[i] == 1 and start is None:
            start = mask_time[i]

        if mask_vals[i] == 0 and start is not None:
            end = mask_time[i]
            intervals.append((start, end))
            start = None

    # close last interval if needed
    if start is not None:
        intervals.append((start, mask_time[-1]))

    return intervals

def _plot_segmented_line_interval(
    ax,
    hybrid_time,
    hybrid_y,
    mask_time,
    mask,
    color,
    label=None,
):
    first_label_used = False
    base_color = color
    for i in range(len(mask)):
        t0 = mask_time[i]
        t1 = mask_time[i + 1] if i < len(mask_time) - 1 else hybrid_time[-1]

        is_dotted = mask[i] == 1
        linestyle = "solid" if is_dotted else "solid"
        color = colors["purple"] if is_dotted else base_color

        # Select hybrid points in [t0, t1)
        idx = np.where((hybrid_time >= t0) & (hybrid_time < t1))[0]

        if len(idx) == 0:
            continue

        lbl = label if not first_label_used else None
        first_label_used = True

        ax.plot(
            hybrid_time[idx],
            hybrid_y[idx],
            color=color,
            linestyle=linestyle,
            label=lbl,
        )

        # Add markers ONLY when switching into or out of loosely dotted
        if is_dotted:
            # Start marker (only if previous interval wasn't loosely dotted)
            if i == 0 or mask[i - 1] == 0:
                ax.scatter(
                    hybrid_time[idx[0]],
                    hybrid_y[idx[0]],
                    color=color,
                    zorder=3,
                    s=7
                )

            # End marker (only if next interval isn't loosely dotted)
            if i == len(mask) - 1 or mask[i + 1] == 0:
                ax.scatter(
                    hybrid_time[idx[-1]],
                    hybrid_y[idx[-1]],
                    color=color,
                    zorder=3,
                    s=7
                )

def _plot_figure(save_path, time_vals, model_used_ts, mean_smm, mean_data, comp_index, comp, num_regions, num_rows, num_cols, figsize, hybrid_label, hybrid_color, plot_smm=True, colname="", ylabel = r"$\mu_I$"):
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    for r in range(num_regions):
        if num_regions == 1:
            ax = axes
        elif num_rows == 1:
            ax = axes[r]
        else:
            ax = axes[r // num_cols, r % num_cols]

        # Extract data
        if(colname != ""):
            col_name = colname
            moment_index = [0 for _ in range(len(compartment_names) * num_regions)]
            moment_index[comp_index + r * len(compartment_names)] = 2
            for m in moment_index:
                col_name += f"{m}"
            y_hybrid = np.sqrt(mean_data[col_name].iloc[:].values)
        else:
            y_hybrid = mean_data.iloc[:, 1 + comp_index + r * len(compartment_names)].values
            
        region_mask = model_used_ts.iloc[:, 1 + r].values

        # Plot SMM
        if plot_smm and mean_smm is not None:
            if(colname != ""):
                col_name = colname
                moment_index = [0 for _ in range(len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y_smm = np.sqrt(mean_smm[col_name].iloc[:].values)
            else:
                y_smm = mean_smm.iloc[:, 1 + comp_index + r * len(compartment_names)].values
            ax.plot(
                mean_smm.Time,
                y_smm,
                color="black",
                linestyle="dashed",
                label="Stochastic",
            )

        # Plot hybrid with segmentation
        mask_time = model_used_ts.iloc[:, 0].values
        region_mask = model_used_ts.iloc[:, 1 + r].values

        _plot_segmented_line_interval(
            ax,
            time_vals,
            y_hybrid,
            mask_time,
            region_mask,
            hybrid_color,
            label=hybrid_label,
        )
        
        # Plot shaded grey background where ODE is used
        mask_time = model_used_ts.iloc[:, 0].values
        mask_vals = model_used_ts.iloc[:, 1 + r].values

        intervals = get_shaded_intervals(mask_time, mask_vals)

        for t0, t1 in intervals:
            ax.axvspan(
                t0,
                t1,
                color="lightgrey",
                alpha=0.5,
                zorder=0,
                linewidth=0,
            )


        # Axis formatting
        if r % num_cols == 0:
            ax.set_ylabel(ylabel)
        # else:
        #     ax.set_yticks([])

        if r // num_cols == num_rows - 1:
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])

        ax.set_title(f"Region {r}")

    # Legend (only once)
    if num_regions == 1:
        ax = axes
    elif num_rows == 1:
        ax = axes[0]
    else:
        ax = axes[0, 0]

    handles, labels = ax.get_legend_handles_labels()

    # Add explanation for loosely dotted segments
    handles.append(plt.Line2D([0], [0], color=colors["purple"], linestyle="solid"))
    labels.append("ODE used")
    handles.append(Patch(facecolor="lightgrey", alpha=0.5, label="ODE used"))
    labels.append("ODE used")

    ax.legend(handles, labels, loc="upper left")

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)


def mean_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, dir_smm):
    comp = list(compartment_colors.keys())[comp_index]

    # Load data
    model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")

    mean_hybrid = pd.read_csv(hybrid_dir + "/joint_mean.csv")
    mean_hybrid_stochastic = pd.read_csv(hybrid_dir + "/stochastic_mean.csv")
    mean_hybrid_deterministic = pd.read_csv(hybrid_dir + "/deterministic_mean.csv")
    time_vals = mean_hybrid.iloc[:, 0].values

    mean_smm = pd.read_csv(dir_smm + "/means.csv") if dir_smm != "" else None

    # Layout calculation
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    # --- Full hybrid ---
    _plot_figure(
        f"{save_dir}/mean_{compartment_names[comp]}_all_regions.png",
        time_vals,
        model_used_ts,
        mean_smm,
        mean_hybrid,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid",
        hybrid_color=compartment_colors[comp][0]
    )

    # --- Stochastic part ---
    _plot_figure(
        f"{save_dir}/mean_{compartment_names[comp]}_all_regions_stoch.png",
        time_vals,
        model_used_ts,
        mean_smm,
        mean_hybrid_stochastic,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid stochastic",
        hybrid_color=compartment_colors[comp][0],
    )

    # --- Deterministic part ---
    _plot_figure(
        f"{save_dir}/mean_{compartment_names[comp]}_all_regions_det.png",
        time_vals,
        model_used_ts,
        mean_smm,
        mean_hybrid_deterministic,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid deterministic",
        hybrid_color=compartment_colors[comp][0],
    )
    
def std_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, dir_smm):
    comp = list(compartment_colors.keys())[comp_index]

    # Load data
    model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")

    moments_hybrid = pd.read_csv(hybrid_dir + "/joint_moments.csv")
    moments_hybrid_stochastic = pd.read_csv(hybrid_dir + "/stochastic_moments.csv")
    moments_hybrid_deterministic = pd.read_csv(hybrid_dir + "/deterministic_moments.csv")
    time_vals = moments_hybrid.iloc[:, 0].values

    moments_smm = pd.read_csv(dir_smm + "/moments.csv") if dir_smm != "" else None

    # Layout calculation
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    # --- Full hybrid ---
    _plot_figure(
        f"{save_dir}/std_{compartment_names[comp]}_all_regions.png",
        time_vals,
        model_used_ts,
        moments_smm,
        moments_hybrid,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid",
        hybrid_color=compartment_colors[comp][0], colname=f"M", ylabel=r"$\sigma_I$"
    )

    # --- Stochastic part ---
    _plot_figure(
        f"{save_dir}/std_{compartment_names[comp]}_all_regions_stoch.png",
        time_vals,
        model_used_ts,
        moments_smm,
        moments_hybrid_stochastic,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid stochastic",
        hybrid_color=compartment_colors[comp][0], colname=f"M", ylabel=r"$\sigma_I$"
    )

    # --- Deterministic part ---
    _plot_figure(
        f"{save_dir}/std_{compartment_names[comp]}_all_regions_det.png",
        time_vals,
        model_used_ts,
        moments_smm,
        moments_hybrid_deterministic,
        comp_index,
        comp,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        hybrid_label="Hybrid deterministic",
        hybrid_color=compartment_colors[comp][0], colname=f"M", ylabel=r"$\sigma_I$"
    )
    

if __name__ == "__main__":
    figsize = (7, 5)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Spatial-Hybrid2/SIR"
    config = "config_4r_10_transm_k2"
    num_regions = 4
    condition = "combined_relation_var_gradient_condition_region"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    tmin = 0
    tmax = 90
    condition_name = ""
    
    smm_dir = f"{dir}/SMM/SIR/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}/"
    save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    # plot_mean_ts_old(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax, tmin, smm_dir, color_smm)
    # plot_var_ts_old(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax, smm_dir, color_smm)
    
    mean_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, smm_dir)
    std_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, smm_dir)
