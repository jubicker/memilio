import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *


def _collect_unique_legend_entries(axes_list):
    """Collect legend handles/labels from multiple axes, deduplicated by label
    (first occurrence wins) and preserving order of first appearance."""
    handles, labels = [], []
    seen = set()
    for ax in axes_list:
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in seen:
                seen.add(li)
                handles.append(hi)
                labels.append(li)
    return handles, labels


def _order_legend_entries_by_condition(handles, labels, base_labels):
    """Reorder legend handles/labels so that, for each base label (e.g. a
    condition name), its deterministic ("... (ODE)") entry is immediately
    followed by its stochastic ("... (stochastic)") entry. Entries not tied
    to any base label keep their original position."""
    label_to_handle = dict(zip(labels, handles))
    emitted = set()
    ordered_labels = []

    for lbl in labels:
        if lbl in emitted:
            continue

        base = next(
            (b for b in base_labels if lbl in (
                b + " (ODE)", b + " (stochastic)")),
            None,
        )

        if base is None:
            ordered_labels.append(lbl)
            emitted.add(lbl)
            continue

        for suffix in (" (ODE)", " (stochastic)"):
            paired = base + suffix
            if paired in label_to_handle and paired not in emitted:
                ordered_labels.append(paired)
                emitted.add(paired)

    ordered_handles = [label_to_handle[l] for l in ordered_labels]
    return ordered_handles, ordered_labels


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


def _plot_figure(save_path, time_vals, model_used_ts, mean_smm, mean_data, comp_index, comp, num_regions, num_rows, num_cols, figsize, hybrid_label, hybrid_color, plot_smm=True, colname="", ylabel=r"$\mu_I$"):
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    for r in range(num_regions):
        if num_regions == 1:
            ax = axes
        elif num_rows == 1:
            ax = axes[r]
        else:
            ax = axes[r // num_cols, r % num_cols]

        # Extract data
        if (colname != ""):
            col_name = colname
            moment_index = [0 for _ in range(
                len(compartment_names) * num_regions)]
            moment_index[comp_index + r * len(compartment_names)] = 2
            for m in moment_index:
                col_name += f"{m}"
            y_hybrid = np.sqrt(mean_data[col_name].iloc[:].values)
        else:
            y_hybrid = mean_data.iloc[:, 1 + comp_index +
                                      r * len(compartment_names)].values

        region_mask = model_used_ts.iloc[:, 1 + r].values

        # Plot SMM
        if plot_smm and mean_smm is not None:
            if (colname != ""):
                col_name = colname
                moment_index = [0 for _ in range(
                    len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y_smm = np.sqrt(mean_smm[col_name].iloc[:].values)
            else:
                y_smm = mean_smm.iloc[:, 1 + comp_index +
                                      r * len(compartment_names)].values
            ax.plot(
                mean_smm.Time,
                y_smm,
                color=colors["dark grey"],
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

        # if r == 1:
        #     ax.set_ylim(bottom=-1000, top=30000)

    # Legend (only once)
    if num_regions == 1:
        ax = axes
    elif num_rows == 1:
        ax = axes[0]
    else:
        ax = axes[num_rows-1, num_cols-1]

    handles, labels = ax.get_legend_handles_labels()

    # Add explanation for loosely dotted segments
    handles.append(plt.Line2D(
        [0], [0], color=colors["orange"], linestyle="solid"))
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
    mean_hybrid_deterministic = pd.read_csv(
        hybrid_dir + "/deterministic_mean.csv")
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
    moments_hybrid_stochastic = pd.read_csv(
        hybrid_dir + "/stochastic_moments.csv")
    moments_hybrid_deterministic = pd.read_csv(
        hybrid_dir + "/deterministic_moments.csv")
    time_vals = moments_hybrid.iloc[:, 0].values

    moments_smm = pd.read_csv(
        dir_smm + "/moments.csv") if dir_smm != "" else None

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
    solid_linestyle="solid",
):
    label_stochastic_used = False
    label_ode_used = False
    for i in range(len(mask)):
        t0 = mask_time[i]
        t1 = mask_time[i + 1] if i < len(mask_time) - 1 else hybrid_time[-1]

        is_dotted = mask[i] == 1
        linestyle = override_linestyle if is_dotted else solid_linestyle
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
    ylabel=r"$\mu_I$",
    plot_smm=True,
    region_titles=None,
):
    # conditions_info: list of dicts with keys:
    #   name, hybrid_dir, ode_color, stoch_color
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    region_axes = []

    filename_smm = "means" if colname == "" else "moments"
    filename_hybrid = "mean" if colname == "" else "moments"

    # Load SMM once if requested
    smm = pd.read_csv(dir_smm + "/" + filename_smm +
                      ".csv") if (dir_smm != "" and plot_smm) else None

    # For each condition load its files
    loaded = []
    for cond in conditions_info:
        hybrid_dir = cond["hybrid_dir"]
        model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")
        hybrid = pd.read_csv(hybrid_dir + "/joint_" + filename_hybrid + ".csv")
        hybrid_stochastic = pd.read_csv(
            hybrid_dir + "/stochastic_" + filename_hybrid + ".csv")
        hybrid_deterministic = pd.read_csv(
            hybrid_dir + "/deterministic_" + filename_hybrid + ".csv")
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
        region_axes.append(ax)

        max = 100000

        # Plot SMM (once)
        if plot_smm and smm is not None:
            if colname != "":
                col_name = colname
                moment_index = [0 for _ in range(
                    len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y_smm = np.sqrt(smm[col_name].iloc[:].values)
            else:
                y_smm = smm.iloc[:, 1 + comp_index +
                                 r * len(compartment_names)].values
            ax.plot(
                smm.Time,
                y_smm,
                color=colors["dark grey"],
                linestyle="dashed",
                label="Stochastic",
            )
            max = y_smm.max() * 1.2

        # For each condition overlay stochastic and deterministic parts
        for cond in loaded:
            model_used_ts = cond["model_used_ts"]
            mask_time = model_used_ts.iloc[:, 0].values
            region_mask = model_used_ts.iloc[:, 1 + r].values

            # stochastic part
            if colname != "":
                col_name = colname
                moment_index = [0 for _ in range(
                    len(compartment_names) * num_regions)]
                moment_index[comp_index + r * len(compartment_names)] = 2
                for m in moment_index:
                    col_name += f"{m}"
                y = np.sqrt(cond["hybrid"][col_name].iloc[:].values)
            else:
                y = cond["hybrid"].iloc[:, 1 + comp_index +
                                        r * len(compartment_names)].values

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
        if region_titles is not None:
            ax.set_title(region_titles[r], pad=11)
        if r % num_cols == 0:
            ax.set_ylabel(ylabel)

        if r // num_cols == num_rows - 1:
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])
        if r == 1:
            ax.set_ylim(bottom=-1000, top=max)

        # if colname == "":
        #     ax.axhline(y=10, color=colors["dark red"],
        #                linestyle="dashed", linewidth=1)

    # Legend (aggregate across all region subplots, since a given condition's
    # ODE/stochastic segment may not appear in every region)
    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)

    handles, labels = _collect_unique_legend_entries(region_axes)
    handles, labels = _order_legend_entries_by_condition(
        handles, labels, [cond["name"] for cond in conditions_info])
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=1)
    fig_leg.tight_layout()
    fig_leg.savefig(save_path + f"_legend.png", dpi=dpi,
                    bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


def mean_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm="", region_titles=None):
    """
    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    region_titles: optional list of titles, one per region, used as subplot titles
    """
    comp = list(compartment_colors.keys())[comp_index]

    # use time from first condition
    first_hybrid = pd.read_csv(
        conditions_info[0]["hybrid_dir"] + "/joint_mean.csv")
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
        region_titles=region_titles,
    )


def std_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm="", region_titles=None):
    """
    region_titles: optional list of titles, one per region, used as subplot titles
    """
    comp = list(compartment_colors.keys())[comp_index]

    first_hybrid = pd.read_csv(
        conditions_info[0]["hybrid_dir"] + "/joint_moments.csv")
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
        region_titles=region_titles,
    )


def _plot_figure_multi_error(
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
    ylabel=r"$|\mu_I-\mu_I^{(stochastic)}|$",
    region_titles=None,
    num_runs=None,
):
    # conditions_info: list of dicts with keys:
    #   name, hybrid_dir, ode_color, stoch_color
    #
    # Plots (on a log scale) the absolute error between each condition's
    # hybrid result and the stochastic (SMM) result. If colname == "" (mean)
    # and num_runs is given, also plots the standard error of the hybrid mean
    # (std/sqrt(num_runs)) as a dark grey line.
    assert dir_smm != "", "dir_smm is required to compute errors against the stochastic model"

    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    region_axes = []

    filename_smm = "means" if colname == "" else "moments"
    filename_hybrid = "mean" if colname == "" else "moments"

    smm = pd.read_csv(dir_smm + "/" + filename_smm + ".csv")
    smm_moments = None
    if colname == "" and num_runs is not None:
        smm_moments = pd.read_csv(dir_smm + "/moments.csv")

    # For each condition load its files
    loaded = []
    for cond in conditions_info:
        hybrid_dir = cond["hybrid_dir"]
        model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")
        hybrid = pd.read_csv(hybrid_dir + "/joint_" + filename_hybrid + ".csv")
        entry = {
            "name": cond["name"],
            "model_used_ts": model_used_ts,
            "hybrid": hybrid,
            "ode_color": cond["ode_color"],
            "stoch_color": cond["stoch_color"],
        }
        loaded.append(entry)

    for r in range(num_regions):
        if num_regions == 1:
            ax = axes
        elif num_rows == 1:
            ax = axes[r]
        else:
            ax = axes[r // num_cols, r % num_cols]
        region_axes.append(ax)

        moment_index = [0 for _ in range(
            len(compartment_names) * num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        moment_col = "M" + "".join(str(m) for m in moment_index)

        if colname != "":
            y_smm = np.sqrt(smm[moment_col].iloc[:].values)
        else:
            y_smm = smm.iloc[:, 1 + comp_index +
                             r * len(compartment_names)].values

        if colname == "" and num_runs is not None:
            std_smm = np.sqrt(smm_moments[moment_col].iloc[:].values)
            time_vals_se = time_vals
            if (len(std_smm) != len(time_vals_se)):
                min_len = min(len(std_smm), len(time_vals_se))
                std_smm = std_smm[:min_len]
                time_vals_se = time_vals[:min_len]
            se = std_smm / np.sqrt(num_runs)
            ax.plot(
                time_vals_se,
                se,
                color=colors["dark grey"],
                linestyle="dotted",
                label="standard error",
            )

        for cond in loaded:
            model_used_ts = cond["model_used_ts"]
            mask_time = model_used_ts.iloc[:, 0].values
            region_mask = model_used_ts.iloc[:, 1 + r].values

            if colname != "":
                y_hybrid = np.sqrt(cond["hybrid"][moment_col].iloc[:].values)
            else:
                y_hybrid = cond["hybrid"].iloc[:, 1 + comp_index +
                                               r * len(compartment_names)].values
            if (len(y_hybrid) != len(y_smm)):
                min_len = min(len(y_hybrid), len(y_smm))
                y_hybrid = y_hybrid[:min_len]
                y_smm = y_smm[:min_len]
            y_error = np.abs(y_hybrid - y_smm)

            if np.isnan(y_error).any():
                continue

            _plot_segmented_line_interval_multi(
                ax,
                time_vals,
                y_error,
                mask_time,
                region_mask,
                cond["ode_color"],
                cond["stoch_color"],
                label=f"{cond['name']}",
                override_linestyle="solid",
            )

        # Axis formatting
        ax.set_yscale("log")
        if region_titles is not None:
            ax.set_title(region_titles[r], pad=11)
        if r % num_cols == 0:
            ax.set_ylabel(ylabel)

        if r // num_cols == num_rows - 1:
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])

    # Legend (aggregate across all region subplots, since a given condition's
    # ODE/stochastic segment may not appear in every region)
    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)

    handles, labels = _collect_unique_legend_entries(region_axes)
    handles, labels = _order_legend_entries_by_condition(
        handles, labels, [cond["name"] for cond in conditions_info])
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center')
    fig_leg.tight_layout()
    fig_leg.savefig(save_path + f"_legend.png", dpi=dpi,
                    bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


def mean_error_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, num_runs, dir_smm="", region_titles=None):
    """
    Plots the absolute error between the hybrid mean and the stochastic (SMM)
    mean (on a log scale), with the standard error of the hybrid mean
    (std/sqrt(num_runs)) plotted as a dark grey line.

    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    region_titles: optional list of titles, one per region, used as subplot titles
    """
    comp = list(compartment_colors.keys())[comp_index]

    # use time from first condition
    first_hybrid = pd.read_csv(
        conditions_info[0]["hybrid_dir"] + "/joint_mean.csv")
    time_vals = first_hybrid.iloc[:, 0].values

    # Layout
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    _plot_figure_multi_error(
        f"{save_dir}/mean_error_{compartment_names[comp]}_multiple_conditions_all_regions",
        time_vals,
        conditions_info,
        comp_index,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        dir_smm,
        colname="",
        ylabel=r"err$(\mu_I)$",
        region_titles=region_titles,
        num_runs=num_runs,
    )
    plt.close("all")


def std_error_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm="", region_titles=None):
    """
    Plots (on a log scale) the absolute error between the hybrid std and the
    stochastic (SMM) std.

    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    region_titles: optional list of titles, one per region, used as subplot titles
    """
    comp = list(compartment_colors.keys())[comp_index]

    first_hybrid = pd.read_csv(
        conditions_info[0]["hybrid_dir"] + "/joint_moments.csv")
    time_vals = first_hybrid.iloc[:, 0].values

    # Layout
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    _plot_figure_multi_error(
        f"{save_dir}/std_error_{compartment_names[comp]}_multiple_conditions_all_regions",
        time_vals,
        conditions_info,
        comp_index,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        dir_smm,
        colname="M",
        ylabel=r"err$(\sigma_I)$",
        region_titles=region_titles,
    )
    plt.close("all")


def mean_error_mean_std_bar_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, dir_smm):
    """
    For each condition in conditions_info, computes the absolute error between
    the hybrid mean and the stochastic (SMM) mean, and the absolute error
    between the hybrid std and the stochastic (SMM) std, aggregated over all
    regions and all time steps, and plots the mean of each error as a grouped
    bar chart (one group per condition, one bar for the mean error and one
    for the std error).

    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    """
    assert dir_smm != "", "dir_smm is required to compute errors against the stochastic model"

    comp = list(compartment_colors.keys())[comp_index]
    n_comp = len(compartment_names)

    smm_mean = pd.read_csv(dir_smm + "/means.csv")
    smm_moments = pd.read_csv(dir_smm + "/moments.csv")

    names = []
    mean_errors = []
    std_errors = []
    bar_colors = []

    for cond in conditions_info:
        hybrid_mean = pd.read_csv(cond["hybrid_dir"] + "/joint_mean.csv")
        hybrid_moments = pd.read_csv(cond["hybrid_dir"] + "/joint_moments.csv")

        mean_errors_r = []
        std_errors_r = []
        for r in range(num_regions):
            y_hybrid = hybrid_mean.iloc[:, 1 + comp_index + r * n_comp].values
            y_smm = smm_mean.iloc[:, 1 + comp_index + r * n_comp].values
            min_len = min(len(y_hybrid), len(y_smm))
            mean_errors_r.append(
                np.abs(y_hybrid[:min_len] - y_smm[:min_len]))

            moment_index = [0 for _ in range(n_comp * num_regions)]
            moment_index[comp_index + r * n_comp] = 2
            moment_col = "M" + "".join(str(m) for m in moment_index)

            y_hybrid_std = np.sqrt(hybrid_moments[moment_col].iloc[:].values)
            y_smm_std = np.sqrt(smm_moments[moment_col].iloc[:].values)
            min_len = min(len(y_hybrid_std), len(y_smm_std))
            std_errors_r.append(
                np.abs(y_hybrid_std[:min_len] - y_smm_std[:min_len]))

        mean_errors_r = np.concatenate(mean_errors_r)
        mean_errors_r = mean_errors_r[~np.isnan(mean_errors_r)]
        std_errors_r = np.concatenate(std_errors_r)
        std_errors_r = std_errors_r[~np.isnan(std_errors_r)]

        names.append(cond["name"])
        mean_errors.append(mean_errors_r.mean())
        std_errors.append(std_errors_r.mean())
        bar_colors.append(cond.get("stoch_color", colors["dark grey"]))

    x = np.arange(len(names))
    width = 0.35

    # Faded fill for the std bars, but keep the hatch lines drawn in the
    # full-opacity bar color (patch alpha would otherwise fade those too).
    std_face_colors = [mcolors.to_rgba(c, alpha=0.5) for c in bar_colors]

    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()

    ax1.bar(x - width / 2, mean_errors, width,
            color=bar_colors, edgecolor=bar_colors)
    ax2.bar(x + width / 2, std_errors, width,
            color=std_face_colors, hatch="//", edgecolor=bar_colors)

    ax1.set_ylabel(r"err$(\mu_I)$")
    ax2.set_ylabel(r"err$(\sigma_I)$")
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)

    legend_handles = [
        Patch(facecolor="grey", edgecolor="grey", label=r"$\mu_I$"),
        Patch(facecolor="grey", alpha=0.5, hatch="//",
              edgecolor="grey", label=r"$\sigma_I$"),
    ]
    ax1.legend(handles=legend_handles, loc="upper right")

    fig.tight_layout()
    fig.savefig(
        f"{save_dir}/mean_error_mean_std_{compartment_names[comp]}_bar_conditions.png", dpi=dpi)
    plt.close(fig)


def _plot_figure_r0_std_mean_ratio_multi(
    save_path,
    time_vals,
    conditions_info,
    comp_index,
    num_regions,
    num_rows,
    num_cols,
    figsize,
    lamdas,
    gamma,
    dir_smm="",
    plot_smm=True,
    region_titles=None,
    r0_ylabel=r"$R_0(\mu)$",
    ratio_ylabel=r"$\sigma_I/\mu_I$",
):
    # conditions_info: list of dicts with keys:
    #   name, hybrid_dir, ode_color, stoch_color
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)

    n_comp = len(compartment_names)

    # Load SMM once if requested
    smm_mean = pd.read_csv(
        dir_smm + "/means.csv") if (dir_smm != "" and plot_smm) else None
    smm_moments = pd.read_csv(
        dir_smm + "/moments.csv") if (dir_smm != "" and plot_smm) else None

    # For each condition load its files
    loaded = []
    for cond in conditions_info:
        hybrid_dir = cond["hybrid_dir"]
        model_used_ts = pd.read_csv(hybrid_dir + "/model_used.csv")
        mean_df = pd.read_csv(hybrid_dir + "/joint_mean.csv")
        moments_df = pd.read_csv(hybrid_dir + "/joint_moments.csv")
        loaded.append({
            "name": cond["name"],
            "model_used_ts": model_used_ts,
            "mean": mean_df,
            "moments": moments_df,
            "ode_color": cond["ode_color"],
            "stoch_color": cond["stoch_color"],
        })

    region_axes = []
    for r in range(num_regions):
        if num_regions == 1:
            ax = axes
        elif num_rows == 1:
            ax = axes[r]
        else:
            ax = axes[r // num_cols, r % num_cols]

        ax2 = ax.twinx()
        region_axes.append((ax, ax2))

        base = 1 + r * n_comp

        moment_index = [0 for _ in range(n_comp * num_regions)]
        moment_index[comp_index + r * n_comp] = 2
        col_name = "M" + "".join(str(m) for m in moment_index)

        # Plot SMM (once)
        if plot_smm and smm_mean is not None and smm_moments is not None:
            S = smm_mean.iloc[:, base + 0].values
            I = smm_mean.iloc[:, base + 1].values
            R = smm_mean.iloc[:, base + 2].values
            y_r0_smm = lamdas[r] * (S) / gamma

            std_smm = np.sqrt(smm_moments[col_name].iloc[:].values)
            mean_smm_comp = smm_mean.iloc[:, 1 +
                                          comp_index + r * n_comp].values
            y_ratio_smm = std_smm / mean_smm_comp

            ax.plot(smm_mean.Time, y_r0_smm, color=colors["dark grey"],
                    linestyle="solid", label=r"Stochastic $R_0(\mu)$")
            ax.set_ylim(bottom=0, top=4.1)
            ax2.plot(smm_moments.Time, y_ratio_smm, color=colors["dark grey"],
                     linestyle="dashed", label=r"Stochastic $\sigma_I/\mu_I$")

        # For each condition overlay R0 (left axis) and std/mean (right axis)
        for cond in loaded:
            model_used_ts = cond["model_used_ts"]
            mask_time = model_used_ts.iloc[:, 0].values
            region_mask = model_used_ts.iloc[:, 1 + r].values

            mean_df = cond["mean"]
            moments_df = cond["moments"]

            S = mean_df.iloc[:, base + 0].values
            y_r0 = lamdas[r] * (S) / gamma

            std = np.sqrt(moments_df[col_name].iloc[:].values)
            mean_comp = mean_df.iloc[:, 1 + comp_index + r * n_comp].values
            y_ratio = std / mean_comp

            _plot_segmented_line_interval_multi(
                ax,
                time_vals,
                y_r0,
                mask_time,
                region_mask,
                cond["ode_color"],
                cond["stoch_color"],
                label=f"{cond['name']} ($R_0(\mu)$)",
                override_linestyle="solid",
                solid_linestyle="solid",
            )

            _plot_segmented_line_interval_multi(
                ax2,
                time_vals,
                y_ratio,
                mask_time,
                region_mask,
                cond["ode_color"],
                cond["stoch_color"],
                label=fr"{cond['name']} ($\sigma_I/\mu_I$)",
                override_linestyle="dashed",
                solid_linestyle="dashed",
            )

        # Reference lines
        ax.axhline(1.0, color=colors['dark red'], linestyle="solid",
                   linewidth=1, label=r"$R_0=1$")
        ax2.axhline(0.7, color=colors['dark red'], linestyle="dashed",
                    linewidth=1, label=r"$\sigma_I/\mu_I=0.7$")
        # ax2.axhline(0.4, color=colors['red'], linestyle="dashed",
        #             linewidth=1, label=r"$\sigma/\mu=0.4$")

        # Axis formatting
        if region_titles is not None:
            ax.set_title(region_titles[r], pad=11)
        if r % num_cols == 0:
            ax.set_ylabel(r0_ylabel)
        if (r + 1) % num_cols == 0 or r == num_regions - 1:
            ax2.set_ylabel(ratio_ylabel)

        if r // num_cols == num_rows - 1:
            ax.set_xlabel("Time [days]")
        else:
            ax.set_xticks([])

    # Legend (aggregate both axes of every subplot, since a given condition's
    # ODE/stochastic segment may not appear in every region)
    all_axes = [ax for pair in region_axes for ax in pair]
    handles, labels = _collect_unique_legend_entries(all_axes)
    base_labels = [
        base
        for cond in loaded
        for base in (f"{cond['name']} ($R_0$)", fr"{cond['name']} ($\sigma/\mu$)")
    ]
    handles, labels = _order_legend_entries_by_condition(
        handles, labels, base_labels)

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)

    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.tight_layout()
    fig_leg.savefig(save_path + f"_legend.png", dpi=dpi,
                    bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


def r0_std_mean_ratio_multiple_conditions(save_dir, conditions_info, comp_index, num_regions, figsize, lamdas, gamma, dir_smm="", region_titles=None):
    """
    conditions_info: list of dicts: {"name":..., "hybrid_dir":..., "ode_color":..., "stoch_color":...}
    region_titles: optional list of titles, one per region, used as subplot titles

    For each region plots, over time:
      - R0(t) = lamdas[r] * (S - I - R) / gamma on the left y-axis
      - the std/mean ratio (sigma/mu) of the `comp_index` compartment on the right y-axis
    Two y-axes are used since both quantities can have very different scales.
    """
    comp = list(compartment_colors.keys())[comp_index]

    if conditions_info:
        # use time from first condition
        first_hybrid = pd.read_csv(
            conditions_info[0]["hybrid_dir"] + "/joint_mean.csv")
        time_vals = first_hybrid.iloc[:, 0].values
    else:
        # no conditions: fall back to the stochastic (SMM) time axis
        assert dir_smm != "", "dir_smm is required when conditions_info is empty"
        time_vals = pd.read_csv(dir_smm + "/means.csv").iloc[:, 0].values

    # Layout
    num_cols = int(np.ceil(np.sqrt(num_regions)))
    num_rows = int(np.ceil(num_regions / num_cols))

    _plot_figure_r0_std_mean_ratio_multi(
        f"{save_dir}/R0_std_mean_ratio_{compartment_names[comp]}_multiple_conditions_all_regions",
        time_vals,
        conditions_info,
        comp_index,
        num_regions,
        num_rows,
        num_cols,
        figsize,
        lamdas,
        gamma,
        dir_smm,
        region_titles=region_titles,
    )


if __name__ == "__main__":
    figsize = (7, 5)
    dir = "/Users/julia/sim_outputs/output"
    save_dir = "/Users/julia/sim_outputs/output"
    hybrid_model = "Spatial-Hybrid2"
    config = "config_SIR_I0_0_1_10_100_no_exchange"
    num_regions = 4
    condition = "combined_relation_var_gradient_condition_region"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    tmin = 0
    condition_name = ""
    # region_titles = [r"Region 1 - $I_{init}=0$, $R_0\approx2$", r"Region 2 - $I_{init}=1$, $R_0\approx2$",
    #                  r"Region 3 - $I_{init}=10$, $R_0\approx2$", r"Region 4 - $I_{init}=100$, $R_0\approx2$"]
    region_titles = [r"$I_{init}^{(1)}=10$, $R_0^{(1)}\approx1$", r"$I_{init}^{(2)}=10$, $R_0^{(2)}\approx1.5$",
                     r"$I_{init}^{(3)}=10$, $R_0^{(3)}\approx2$", r"$I_{init}^{(4)}=10$, $R_0^{(4)}\approx4$"]
    gamma = 1/7.
    lamdas = [0.000001432, 0.00000215, 0.00000286, 0.00000572]
    num_runs = 10000

    smm_dir = f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    # save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    save_dir = f"{save_dir}/{hybrid_model}/{config}"
    os.makedirs(save_dir, exist_ok=True)

    comp_index = 1

    # mean_all_regions(save_dir, hybrid_dir, comp_index,
    #                  num_regions, figsize, smm_dir)
    # std_all_regions(save_dir, hybrid_dir, comp_index,
    #                 num_regions, figsize, smm_dir)

    conditions1 = [
        {"name": "Deterministic", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/pure_ode/{closure_method}/{closure_order}",
         "ode_color": colors['middle green'], "stoch_color": colors['middle green']},
        {"name": r"$\tau_{\mu_I}$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/abs_threshold_region/{closure_method}/{closure_order}",
            "ode_color": colors['rose'], "stoch_color": colors['purple']},
        {"name": "combined", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/combined_region/{closure_method}/{closure_order}",
            "ode_color": colors['middle blue'], "stoch_color": colors['dark blue']},
        {"name": r"$\sigma_I/\mu_I$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/mean_stddev_relation_region/{closure_method}/{closure_order}",
            "ode_color": colors['orange'], "stoch_color": colors['brown']},
        {"name": r"$R_0(\mu_I)$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/R0_region/{closure_method}/{closure_order}",
            "ode_color": colors['light teal'], "stoch_color": colors['teal']},
    ]

    conditions2 = [
        {"name": "pure_ode", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/pure_ode/{closure_method}/{closure_order}",
            "ode_color": colors['orange'], "stoch_color": colors['brown']},
        {"name": "pure_stochastic", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/pure_stochastic/{closure_method}/{closure_order}",
            "ode_color": colors['rose'], "stoch_color": colors['purple']}
    ]

    # mean_multiple_conditions(
    #     save_dir, conditions2, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)
    # std_multiple_conditions(
    #     save_dir, conditions2, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)
    # mean_error_multiple_conditions(
    #     save_dir, conditions2, comp_index, num_regions, figsize, num_runs, smm_dir, region_titles=region_titles)
    # std_error_multiple_conditions(
    #     save_dir, conditions2, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)

    # mean_multiple_conditions(
    #     save_dir, conditions1, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)
    # std_multiple_conditions(
    #     save_dir, conditions1, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)
    # mean_error_multiple_conditions(
    #     save_dir, conditions1, comp_index, num_regions, figsize, num_runs, smm_dir, region_titles=region_titles)
    # std_error_multiple_conditions(
    #     save_dir, conditions1, comp_index, num_regions, figsize, smm_dir, region_titles=region_titles)
    mean_error_mean_std_bar_conditions(
        save_dir, conditions1, comp_index, num_regions, figsize, smm_dir)

    # r0_std_mean_ratio_multiple_conditions(
    #     save_dir, [], comp_index, num_regions, figsize, lamdas, gamma, smm_dir, region_titles=region_titles)
