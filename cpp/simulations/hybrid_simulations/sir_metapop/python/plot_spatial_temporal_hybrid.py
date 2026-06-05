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
        color = colors["orange"] if is_dotted else base_color

        # Select hybrid points in [t0, t1)
        idx = np.where((hybrid_time >= t0) & (hybrid_time < t1))[0]

        if len(idx) == 0:
            continue
        
        if color != base_color:
            label = "Hybrid (ODE)"
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
    R0_values = [r"$R_0\approx 4$", r"$R_0\approx 4$", r"$R_0\approx 4$", r"$R_0\approx 4$"]
    I0_values = [r"$I_0=0$", r"$I_0=1$", r"$I_0=10$", r"$I_0=100$"]
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
                alpha=0.0,
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

        ax.set_title(f"Region {r}, " + I0_values[r] + ", " + R0_values[r])
        if r==1:
            ax.set_ylim(bottom = -1000, top = 30000)

    # Legend (only once)
    if num_regions == 1:
        ax = axes
    elif num_rows == 1:
        ax = axes[0]
    else:
        ax = axes[num_rows-1, num_cols-1]

    handles, labels = ax.get_legend_handles_labels()

    # Add explanation for loosely dotted segments
    handles.append(plt.Line2D([0], [0], color=colors["orange"], linestyle="solid"))
    labels.append("Hybrid (ODE)")
    # handles.append(Patch(facecolor="lightgrey", alpha=0.5, label="ODE used"))
    # labels.append("ODE used")

    ax.legend(handles, labels, loc="upper right")

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
        hybrid_label="Hybrid (stochastic)",
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
        hybrid_label="Hybrid (stochastic)",
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
        hybrid_label="Hybrid(deterministic)",
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
        hybrid_label="Hybrid (stochastic)",
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
        hybrid_label="Hybrid (stochastic)",
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
        hybrid_label="Hybrid (stochastic)",
        hybrid_color=compartment_colors[comp][0], colname=f"M", ylabel=r"$\sigma_I$"
    )

def _plot_segmented_line_interval_multi(
    ax,
    hybrid_time,
    hybrid_y,
    mask_time,
    mask,
    color_dotted,
    color_solid,
    label,
    override_linestyle,
):
    label_stochastic_used = False
    label_ode_used = False
    for i in range(len(mask)):
        t0 = mask_time[i]
        t1 = mask_time[i + 1] if i < len(mask_time) - 1 else hybrid_time[-1]

        is_dotted = mask[i] == 1
        linestyle = override_linestyle if is_dotted else "solid"
        lbl = None
        if is_dotted and not label_ode_used:
            lbl = label + " (ODE)"
            label_ode_used = True
        elif not is_dotted and not label_stochastic_used:
            lbl = label + " (stochastic)"
            label_stochastic_used = True

        # Select hybrid points in [t0, t1)
        idx = np.where((hybrid_time >= t0) & (hybrid_time < t1))[0]

        if len(idx) == 0:
            continue
        if is_dotted:
            ax.plot(
                hybrid_time[idx],
                hybrid_y[idx],
                color=color_dotted,
                linestyle=linestyle,
                label=lbl
            )
        else:
            ax.plot(
                hybrid_time[idx],
                hybrid_y[idx],
                color=color_solid,
                linestyle=linestyle,
                label=lbl
            )

        # Add markers ONLY when switching into or out of loosely dotted
        if is_dotted:
            # Start marker (only if previous interval wasn't loosely dotted)
            if i == 0 or mask[i - 1] == 0:
                ax.scatter(
                    hybrid_time[idx[0]],
                    hybrid_y[idx[0]],
                    color=color_dotted,
                    zorder=3,
                    s=7
                )

            # End marker (only if next interval isn't loosely dotted)
            if i == len(mask) - 1 or mask[i + 1] == 0:
                ax.scatter(
                    hybrid_time[idx[-1]],
                    hybrid_y[idx[-1]],
                    color=color_dotted,
                    zorder=3,
                    s=7
                )

def _plot_figure_multi(
    save_path,
    time_vals,
    conditions_info,
    comp_index,
    num_regions,
    num_rows,
    num_cols,
    figsize,
    dir_smm,
    colname="",
    ylabel = r"$\mu_I$",
    plot_smm=True,
):
    # conditions_info: list of dicts with keys:
    #   name, hybrid_dir, ode_color, stoch_color
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    
    filename_smm = "means" if colname == "" else "moments"
    filename_hybrid = "mean" if colname == "" else "moments"

    # Load SMM once if requested
    smm = pd.read_csv(dir_smm + "/" + filename_smm + ".csv") if (dir_smm != "" and plot_smm) else None

    # For each condition load its files
    loaded = []
    for cond in conditions_info:
        hybrid_dir = cond["hybrid_dir"]
        model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")
        hybrid = pd.read_csv(hybrid_dir + "/joint_" + filename_hybrid + ".csv")
        hybrid_stochastic = pd.read_csv(hybrid_dir + "/stochastic_" + filename_hybrid + ".csv")
        hybrid_deterministic = pd.read_csv(hybrid_dir + "/deterministic_" + filename_hybrid + ".csv")
        loaded.append({
            "name": cond["name"],
            "model_used_ts": model_used_ts,
            "hybrid": hybrid,
            "stoch": hybrid_stochastic,
            "det": hybrid_deterministic,
            "ode_color": cond["ode_color"],
            "stoch_color": cond["stoch_color"],
        })

    for r in range(num_regions):
        if num_regions == 1:
            ax = axes
        elif num_rows == 1:
            ax = axes[r]
        else:
            ax = axes[r // num_cols, r % num_cols]
            
        max = 100000

        # Plot SMM (once)
        if plot_smm and smm is not None:
            if colname != "":
                col_name = colname
                moment_index = [0 for _ in range(len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y_smm = np.sqrt(smm[col_name].iloc[:].values)
            else:
                y_smm = smm.iloc[:, 1 + comp_index + r * len(compartment_names)].values
            ax.plot(
                smm.Time,
                y_smm,
                color="black",
                linestyle="dashed",
                label="Stochastic",
            )
            max = y_smm.max() * 1.5

        # For each condition overlay stochastic and deterministic parts
        for cond in loaded:
            model_used_ts = cond["model_used_ts"]
            mask_time = model_used_ts.iloc[:, 0].values
            region_mask = model_used_ts.iloc[:, 1 + r].values

            # stochastic part
            if colname != "":
                col_name = colname
                moment_index = [0 for _ in range(len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y = np.sqrt(cond["hybrid"][col_name].iloc[:].values)
            else:
                y = cond["hybrid"].iloc[:, 1 + comp_index + r * len(compartment_names)].values

            # stochastic (solid)
            _plot_segmented_line_interval_multi(
                ax,
                time_vals,
                y,
                mask_time,
                region_mask,
                cond["ode_color"],
                cond["stoch_color"],
                label=f"{cond['name']}",
                override_linestyle="solid",
            )

        # Axis formatting
        if r % num_cols == 0:
            ax.set_ylabel(ylabel)

        if r // num_cols == num_rows - 1:
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])
        if r==1:
            ax.set_ylim(bottom = -1000, top = max)

        I0_values = [r"$I_0=0$", r"$I_0=1$", r"$I_0=10$", r"$I_0=100$"]
        R0_values = [r"$R_0\approx 4$"] * 4
        ax.set_title(f"Region {r}, " + I0_values[r] + ", " + R0_values[r])

    # Legend
    if num_regions == 1:
        ax = axes
    elif num_rows == 1:
        ax = axes[0]
    else:
        ax = axes[num_rows-1, num_cols-1]

    #ax.legend(handles, labels, loc="upper right")

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)
    
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center')      
    fig_leg.tight_layout()
    fig_leg.savefig(save_path + f"_legend.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


def mean_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm=""):
    """
    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    """
    comp = list(compartment_colors.keys())[comp_index]

    # use time from first condition
    first_hybrid = pd.read_csv(conditions_info[0]["hybrid_dir"] + "/joint_mean.csv")
    time_vals = first_hybrid.iloc[:, 0].values

    # Layout
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    _plot_figure_multi(
        f"{save_dir}/mean_{compartment_names[comp]}_multiple_conditions_all_regions",
        time_vals,
        conditions_info,
        comp_index,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        dir_smm,
        colname="",
        ylabel=r"$\mu_I$",
    )


def std_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm=""):
    comp = list(compartment_colors.keys())[comp_index]

    first_hybrid = pd.read_csv(conditions_info[0]["hybrid_dir"] + "/joint_moments.csv")
    time_vals = first_hybrid.iloc[:, 0].values

    # Layout
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    _plot_figure_multi(
        f"{save_dir}/std_{compartment_names[comp]}_multiple_conditions_all_regions",
        time_vals,
        conditions_info,
        comp_index,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        dir_smm,
        colname="M",
        ylabel=r"$\sigma_I$",
    )

 
if __name__ == "__main__":
    figsize = (7, 5)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Spatial-Hybrid2"
    config = "config_SIR_I0_0_1_10_100_no_exchange"
    num_regions = 4
    condition = "fixed_tp_region"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    tmin = 0
    condition_name = ""
    
    smm_dir = f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    
    # mean_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, smm_dir)
    # std_all_regions(save_dir, hybrid_dir, comp_index, num_regions, figsize, smm_dir)
    
    conditions_info = [
        {"name": "pure_ode", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/pure_ode/{closure_method}/{closure_order}", "ode_color": colors['orange'], "stoch_color": colors['brown']},
        {"name": "fixed_tp", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/fixed_tp_region/{closure_method}/{closure_order}", "ode_color": colors['rose'], "stoch_color": colors['purple']},
        {"name": "abs_threshold", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/abs_threshold_condition_region/{closure_method}/{closure_order}", "ode_color": colors['middle blue'], "stoch_color": colors['dark blue']},
        {"name": "relation_var_gradient", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/combined_relation_var_gradient_condition_region/{closure_method}/{closure_order}", "ode_color": colors['middle green'], "stoch_color": colors['dark green']}
    ]
    
    mean_multiple_conditions(f"H:/Documents/TemporalHybridModel/Spatial-Hybrid2/{config}", conditions_info, comp_index, num_regions, figsize, smm_dir)
    std_multiple_conditions(f"H:/Documents/TemporalHybridModel/Spatial-Hybrid2/{config}", conditions_info, comp_index, num_regions, figsize, smm_dir)
