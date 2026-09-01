import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from settings import *
import os


def plot_scenarios_legend(scenarios, base_color, scenario_start_date, base_dir, figsize, id):
    handles = [
        Line2D([0], [0], marker='x', color='black', linestyle='None',
               label='Real incidence'),
    ]
    handles.append(
        Line2D([0], [0], color=colors["dark grey"], label=r'$\mu_I$'))
    handles.append(
        Patch(facecolor=colors["dark grey"], alpha=0.5,
              label=r'$\mu_I \pm 2\sigma_I$'))
    handles.append(Patch(facecolor=base_color,
                         label=f"Simulated data\nuntil {scenario_start_date.date()}"))
    for sc in range(len(scenarios["Path"])):
        handles.append(
            Patch(facecolor=scenarios["Color"][sc],
                  label=scenarios["Name"][sc]))

    fig_legend = plt.figure(figsize=(figsize[0]*1.5, figsize[1]))
    fig_legend.legend(handles=handles, loc='center', ncol=4, frameon=False)
    fig_legend.savefig(base_dir + f"legend_scenarios_{id}.png",
                       bbox_inches='tight')
    plt.close(fig_legend)


def plot_scenarios(
        real_data_file,
        base_dir,
        scenarios,
        base_color,
        num_days,
        scenario_start,
        num_regions,
        region,
        region_index, start_date, pop_size, figsize, intervention_tps=None):

    var_col = ["0" for _ in range(num_regions * 3)]
    var_col[3 * region_index + 1] = "2"
    var_col = "M" + "".join(var_col)
    scenario_start_date = start_date + pd.Timedelta(days=scenario_start)
    scenario_end_date = start_date + pd.Timedelta(days=num_days)

    fig, ax = plt.subplots(figsize=figsize)
    fig2, ax2 = plt.subplots(figsize=(0.7*figsize[0], figsize[1]))

    for sc in range(len(scenarios["Path"])):
        # Read means and moments
        mean_df = pd.read_csv(base_dir + scenarios["Path"][sc] + "/means.csv")
        moment_df = pd.read_csv(
            base_dir + scenarios["Path"][sc] + "/moments.csv")

        mean_df["Date"] = start_date + \
            pd.to_timedelta(mean_df.Time, unit="D")
        moment_df["Date"] = start_date + \
            pd.to_timedelta(moment_df.Time, unit="D")

        time_until_scenario = mean_df[mean_df.Date < scenario_start_date].Date
        time_from_scenario = mean_df[mean_df.Date >= scenario_start_date].Date
        mean_until_scenario = mean_df[mean_df.Date <
                                      scenario_start_date].iloc[:, 1 + region_index*3+1]
        mean_from_scenario = mean_df[mean_df.Date >=
                                     scenario_start_date].iloc[:, 1 + region_index*3+1]
        std_until_scenario = np.sqrt(
            moment_df[moment_df.Date < scenario_start_date][var_col].iloc[:])
        std_from_scenario = np.sqrt(
            moment_df[moment_df.Date >= scenario_start_date][var_col].iloc[:])

        # Plot outputs until scenario start
        ax.fill_between(time_until_scenario, np.maximum(mean_until_scenario - 2*std_until_scenario, 0), mean_until_scenario + 2*std_until_scenario,
                        color=base_color, alpha=0.2, zorder=1)
        ax.plot(time_until_scenario, mean_until_scenario,
                color=base_color, zorder=2)

        # Plot outputs from scenario start until end
        ax.fill_between(time_from_scenario, np.maximum(mean_from_scenario - 2*std_from_scenario, 0), mean_from_scenario + 2*std_from_scenario,
                        color=scenarios["Color"][sc], alpha=0.4, zorder=1)
        ax.plot(time_from_scenario, mean_from_scenario,
                color=scenarios["Color"][sc], zorder=2)

        # Same "from scenario start" outputs on the zoomed-in figure
        ax2.fill_between(time_from_scenario, np.maximum(mean_from_scenario - 2*std_from_scenario, 0), mean_from_scenario + 2*std_from_scenario,
                         color=scenarios["Color"][sc], alpha=0.4, zorder=1)
        ax2.plot(time_from_scenario, mean_from_scenario,
                 color=scenarios["Color"][sc], zorder=2)

    real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])

    ax1 = ax.twinx()
    filtered_df = real_df[(real_df["Datum"] >= start_date)
                          & (real_df["Datum"] < scenario_end_date) & (real_df["Region"] == region) & (real_df["Altersgruppe"] == "00+")]
    ax1.scatter(filtered_df["Datum"], filtered_df["Inzidenz"]/100000*pop_size,
                marker='x', color='black', label='real incidence')
    ax1.set_ylabel("Real incidence")

    y_min = min(ax.get_ylim()[0], ax1.get_ylim()[0])
    y_max = max(ax.get_ylim()[1], ax1.get_ylim()[1])
    ax.set_ylim(y_min, y_max)
    ax1.set_ylim(y_min, y_max)

    ax.set_xlabel("Date")
    ax.set_ylabel("Infected [#]")
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    ax1.grid()
    fig.tight_layout()
    fig.savefig(base_dir + f"mean_std_{region}.png")
    plt.close(fig)

    # Zoomed-in figure: same plot, but only from scenario_start onwards
    ax3 = ax2.twinx()
    filtered_df_from_scenario = real_df[(real_df["Datum"] >= scenario_start_date)
                                        & (real_df["Datum"] < scenario_end_date) & (real_df["Region"] == region) & (real_df["Altersgruppe"] == "00+")]
    ax3.scatter(filtered_df_from_scenario["Datum"], filtered_df_from_scenario["Inzidenz"]/100000*pop_size,
                marker='x', color='black', label='real incidence')
    ax3.set_ylabel("Real incidence")

    y_min2 = min(ax2.get_ylim()[0], ax3.get_ylim()[0])
    y_max2 = max(ax2.get_ylim()[1], ax3.get_ylim()[1])
    ax2.set_ylim(y_min2, y_max2)
    ax3.set_ylim(y_min2, y_max2)

    ax2.set_xlabel("Date")
    ax2.set_ylabel("Infected [#]")
    ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")

    if intervention_tps:
        for tp in intervention_tps:
            tp_date = start_date + pd.Timedelta(days=tp)
            if scenario_start_date <= tp_date < scenario_end_date:
                ax2.axvline(tp_date, color=colors["dark grey"],
                            linestyle="--", zorder=3)

    ax3.grid()
    fig2.tight_layout()
    fig2.savefig(base_dir + f"mean_std_{region}_from_scenario.png")
    plt.close(fig2)

    plot_scenarios_legend(scenarios, base_color,
                          scenario_start_date, base_dir, figsize, id=region)


def plot_scenario_names_legend(scenarios, base_dir, figsize, id):
    handles = [
        Patch(facecolor=scenarios["Color"][sc], label=scenarios["Name"][sc])
        for sc in range(len(scenarios["Path"]))]

    fig_legend = plt.figure(figsize=figsize)
    fig_legend.legend(handles=handles, loc='center', ncol=2, frameon=False)
    fig_legend.savefig(base_dir + f"legend_scenario_names_{id}.png",
                       bbox_inches='tight')
    plt.close(fig_legend)


def plot_model_usage_and_runtime(base_dir, scenarios, figsize, id, smm_runtime=0.):
    percentages = []
    runtimes = []

    for sc in range(len(scenarios["Path"])):
        model_used_df = pd.read_csv(
            base_dir + scenarios["Path"][sc] + "/model_used.csv")
        # First column is time, the remaining columns hold 0/1 per region
        model_used_values = model_used_df.iloc[:, 1:]
        percentages.append(
            (model_used_values.size - model_used_values.values.sum()) / model_used_values.size * 100)

        runtime_df = pd.read_csv(
            base_dir + scenarios["Path"][sc] + "/total_time.csv")
        runtimes.append(runtime_df["Runtime"].iloc[0])

    labels = list(scenarios["Name"])
    bar_colors = list(scenarios["Color"])

    has_smm = smm_runtime != 0
    x = np.arange(len(labels))
    width = 0.35

    if has_smm:

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0.15,
                           1, 0.75])
        ax_smm = fig.add_axes([1.1, 0.15,
                               0.25, 0.75])
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax_smm = None
    ax2 = ax.twinx()

    ax.bar(x - width/2, percentages, width, color=bar_colors,
           edgecolor=bar_colors)
    runtime_bars = ax2.bar(x + width/2, runtimes, width, hatch="//")
    for patch, c in zip(runtime_bars, bar_colors):
        patch.set_facecolor(to_rgba(c, alpha=0.2))
        patch.set_edgecolor(to_rgba(c, alpha=1.0))

    ax.set_ylabel("Stochastic model usage [%]")
    ax.set_ylim(0, 100)
    if (ax_smm is None):
        ax2.set_ylabel("Runtime [s]")
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])

    if ax_smm is not None:
        x_smm = np.array([0])
        width_smm = 0.1
        ax_smm.bar(x_smm - width_smm/2, [100], width_smm,
                   color=colors["middle grey"], edgecolor=colors["middle grey"])
        ax2_smm = ax_smm.twinx()
        ax2_smm.bar(x_smm + width_smm/2, [smm_runtime], width_smm,
                    color=to_rgba(colors["middle grey"], alpha=0.2),
                    edgecolor=to_rgba(colors["middle grey"], alpha=1.0), hatch="//")
        ax_smm.set_ylim(0, 100)
        ax_smm.set_yticklabels([])
        ax2_smm.set_ylabel("Runtime [s]")
        ax_smm.set_xticks(x_smm)
        ax_smm.set_xticklabels([""])

    handles = [Patch(facecolor=colors["dark grey"], edgecolor=colors["dark grey"],
                     label="Stochastic model usage [%]"),
               Patch(facecolor=to_rgba(colors["dark grey"], alpha=0.2),
                     edgecolor=to_rgba(colors["dark grey"], alpha=1.0),
                     hatch="//", label="Runtime [s]")]
    ax.legend(handles=handles, loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)

    if ax_smm is None:
        fig.tight_layout()
    fig.savefig(base_dir + f"model_usage_and_runtime_{id}.png",
                dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    if has_smm:
        labels = labels + ["Base\n(fully stochastic)"]
        bar_colors = bar_colors + [colors["middle grey"]]

    plot_scenario_names_legend(
        {"Path": labels, "Name": labels, "Color": bar_colors}, base_dir, figsize, id=id)


if __name__ == "__main__":
    real_data_df = "/Users/julia/repos/grippeweb_data/ILI_df_regions.csv"
    real_data_df_germany = "/Users/julia/repos/grippeweb_data/ILI_df_germany.csv"
    base_dir = "/Users/julia/sim_outputs/100ppa_jureca/"
    regions = ["Norden (West)", "Osten", "Sueden", "Mitte (West)"]
    regions_germany = ["Bundesweit"]
    base_color = colors['middle blue']
    scenario_start = 1323
    num_days = 5 * 365 - 4
    start_date = pd.to_datetime("2016-08-01")
    pop_sizes = [int(135314.52), int(160939.05),
                 int(244948.26), int(294569.57)]
    pop_size_germany = [int(835771.40)]
    intervention_tps = [1323, 1353, 1554, 1577, 1680, 1736]
    condition = "combined_region/truncation/closure_order_3"
    scenarios = {
        "Path": [
            f"config_influenza_regions_full_base_32c/{condition}",
            f"config_influenza_regions_full_global_32c/{condition}",
            f"config_influenza_regions_full_regional_no_32c/{condition}",
            f"config_influenza_regions_full_regional_strict_32c/{condition}",
            f"config_influenza_regions_full_total_lifting_32c/{condition}"],
        "Name": [
            "Base",
            "Global NPIs",
            "No NPIs\nin the South",
            "Stricter\nNPIs in the South",
            "Global NPIs\nwith lifting"],
        "Color": [
            colors["dark blue"],
            colors["purple"],
            colors["orange"],
            colors["teal"],
            colors["dark green"]]}
    scenarios_germany = {
        "Path": [
            f"config_influenza_germany_full_base_32c/{condition}",
            f"config_influenza_germany_full_global_32c/{condition}",
            f"config_influenza_germany_full_total_lifting_32c/{condition}"
        ],
        "Name": [
            "Base",
            "Global NPIs",
            "Global NPIs\nwith lifting"
        ],
        "Color": [
            colors["dark teal"],
            colors["rose"],
            colors["middle green"],
        ]}

    figsize = (9, 4)
    fig_size_bar = (0.5*figsize[0], 1.5*figsize[1])

    # for r, region in enumerate(regions):
    #     plot_scenarios(real_data_df, base_dir, scenarios, base_color, num_days, scenario_start, num_regions=len(
    #         regions), region=region, region_index=r, start_date=start_date, pop_size=pop_sizes[r], figsize=figsize, intervention_tps=intervention_tps)

    # plot_model_usage_and_runtime(
    #     base_dir, scenarios, fig_size_bar, id="regions")

    # for r, region in enumerate(regions_germany):
    #     plot_scenarios(real_data_df_germany, base_dir, scenarios_germany, base_color, num_days, scenario_start, num_regions=len(
    #         regions_germany), region=region, region_index=r, start_date=start_date, pop_size=pop_size_germany[r], figsize=figsize, intervention_tps=intervention_tps)

    plot_model_usage_and_runtime(
        base_dir, scenarios_germany, fig_size_bar, id="germany", smm_runtime=186)
