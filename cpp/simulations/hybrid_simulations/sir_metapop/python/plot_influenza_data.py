import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from settings import *


def plot_new_infections(sim_output_filepath, save_dir, real_data_file, start_date, tmax):
    sim_outputs = {"Date": [], "Mean": []}
    df = pd.read_csv(sim_output_filepath + "/new_infections_mean.csv")
    current_date = start_date
    for week in range(0, len(df), 7):
        incidence = df.C1.iloc[week: week + 7].sum()
        sim_outputs["Date"].append(current_date)
        sim_outputs["Mean"].append(incidence)
        current_date += pd.Timedelta(days=7)
    sim_outputs = pd.DataFrame(sim_outputs)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            marker='o', linestyle='-', color='blue', label='simulated')

    if (real_data_file != ""):
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        filtered_df = real_df[(real_df["Datum"] >= start_date)
                              & (real_df["Datum"] <= end_date)]
        ax.scatter(filtered_df["Datum"], filtered_df["Inzidenz"],
                   marker='x', color='black', label='real')

    ax.set_xlabel("Date")
    ax.set_ylabel("Incidence")
    ax.set_xticks(sim_outputs["Date"].iloc[::8])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, loc='upper left')
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir +
                f"incidence.png")


def plot_infected(sim_output_filepath, save_dir, start_date, p_lower, p_upper):
    sim_outputs = {"Date": [], "Mean": [], "Std": [], "Lower": [], "Upper": []}
    df_mean = pd.read_csv(sim_output_filepath + "/means.csv")
    df_moments = pd.read_csv(sim_output_filepath + "/moments.csv")
    df_lower = pd.read_csv(sim_output_filepath + f"/{p_lower}.csv")
    df_upper = pd.read_csv(sim_output_filepath + f"/{p_upper}.csv")
    current_date = start_date

    for t in range(0, len(df_mean)):
        sim_outputs["Date"].append(current_date + pd.to_timedelta(t, unit="D"))
        sim_outputs["Mean"].append(df_mean.C2.iloc[t])
        sim_outputs["Std"].append(np.sqrt(df_moments.M020.iloc[t]))
        sim_outputs["Lower"].append(df_lower.C2.iloc[t])
        sim_outputs["Upper"].append(df_upper.C2.iloc[t])

    sim_outputs = pd.DataFrame(sim_outputs)
    interval = int(8*7/0.1)

    # Plot mean +- std
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            color=colors['dark red'], label=r'$\mu$')
    ax.fill_between(sim_outputs["Date"], sim_outputs["Mean"] - sim_outputs["Std"], sim_outputs["Mean"] + sim_outputs["Std"],
                    color=colors['dark red'], alpha=0.5, label=r'$\mu \pm \sigma$')
    ax.set_xlabel("Date")
    ax.set_ylabel("Infected [#]")
    ax.set_xticks(sim_outputs["Date"].iloc[::interval])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, loc='upper left')
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir + f"mean_std.png")

    # Plot mean and percentiles
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            color=colors['dark red'], label=r'$\mu$')
    ax.fill_between(sim_outputs["Date"], sim_outputs["Lower"], sim_outputs["Upper"],
                    color=colors['dark red'], alpha=0.5, label=r'$90$-quantile')
    ax.set_xlabel("Date")
    ax.set_ylabel("Infected [#]")
    ax.set_xticks(sim_outputs["Date"].iloc[::interval])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, loc='upper left')
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir + f"percentiles.png")


plot_new_infections("V:/bick_ju/TemporalHybrid/SMM/config_influenza_germany", "",
                    "C:/Users/bick_ju/Documents/repos/grippeweb_data/ILI_df_Germany.csv", pd.to_datetime("2016-08-01"), 3 * 365)
plot_infected("V:/bick_ju/TemporalHybrid/SMM/config_influenza_germany",
              "", pd.to_datetime("2016-08-01"), "p05", "p95")
