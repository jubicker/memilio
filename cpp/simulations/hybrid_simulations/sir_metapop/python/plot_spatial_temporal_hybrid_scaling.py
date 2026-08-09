import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import glob
import numpy as np
from matplotlib.patches import Patch
from settings import *


def parse_scaling_out(filepath):
    """Parse a *.out file with blocks of the form
    'Run with <N> cores.' followed (a few lines later) by a line
    containing '... time: <seconds>'. Returns a dict {num_cores: time}."""
    times = {}
    current_cores = None
    with open(filepath) as f:
        for line in f:
            m_cores = re.search(r"Run with (\d+) cores?", line)
            if m_cores:
                current_cores = int(m_cores.group(1))
                continue
            m_time = re.search(
                r"time:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line, re.IGNORECASE)
            if m_time and current_cores is not None:
                times[current_cores] = float(m_time.group(1))
                current_cores = None
    return times


def load_scaling_times(dir_path):
    """Find the (single) *.out file in dir_path and parse its timings."""
    out_files = glob.glob(os.path.join(dir_path, "*.out"))
    if not out_files:
        raise FileNotFoundError(f"No .out file found in {dir_path}")
    return parse_scaling_out(out_files[0])


def plot_core_scaling(save_path, models, figsize):
    """models: list of dicts {"name":..., "dir":..., "color":...}"""
    fig, ax = plt.subplots(figsize=figsize)

    all_cores = set()
    for model in models:
        times = load_scaling_times(model["dir"])
        cores = sorted(times.keys())
        all_cores.update(cores)
        runtimes = [times[c] for c in cores]
        ax.plot(cores, runtimes, label=model["name"],
                color=model["color"], marker="o")

    ax.set_yscale("log")
    ax.set_xscale("log", base=2)
    x_ticks = sorted(all_cores)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(c) for c in x_ticks])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Cores [#]")
    ax.set_ylabel("Runtime [s]")
    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)

    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=1)
    fig_leg.tight_layout()
    fig_leg.savefig(save_path + "_legend.png", dpi=dpi,
                    bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


if __name__ == "__main__":
    figsize = (7, 5)
    dir = "/Users/julia/sim_outputs/output"
    save_dir = "/Users/julia/sim_outputs/output"
    hybrid_model = "Spatial-Hybrid2"
    config = "config_SIR_I0_0_1_10_100_no_exchange"
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']

    smm_dir = f"{dir}/SMM/{config}"
    save_dir = f"{save_dir}/{hybrid_model}/{config}"
    os.makedirs(save_dir, exist_ok=True)

    conditions1 = [
        {"name": r"$\tau_{\mu_I}$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/abs_threshold_region/{closure_method}/{closure_order}",
            "color": colors['purple']},
        {"name": "combined", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/combined_region/{closure_method}/{closure_order}",
            "color": colors['dark blue']},
        {"name": r"$\sigma_I/\mu_I$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/mean_stddev_relation_region/{closure_method}/{closure_order}",
            "color": colors['brown']},
        {"name": r"$R_0(\mu_I)$", "hybrid_dir": f"{dir}/{hybrid_model}/{config}/R0_region/{closure_method}/{closure_order}",
            "color": colors['teal']},
        {"name": "Deterministic", "hybrid_dir": f"{dir}/Moments/{config}/{closure_method}/{closure_order}/0.000000",
         "color": colors['middle green']}
    ]

    models = [{"name": "Stochastic", "dir": smm_dir, "color": color_smm}]
    models += [{"name": cond["name"], "dir": cond["hybrid_dir"],
                "color": cond["color"]} for cond in conditions1]

    plot_core_scaling(f"{save_dir}/core_scaling.png", models, figsize)
