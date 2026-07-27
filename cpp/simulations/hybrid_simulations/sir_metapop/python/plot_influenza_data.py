import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from settings import *
import os


def plot_new_infections(sim_output_filepath, save_dir, real_data_file, start_date, tmax, region, pop_size, region_name):
    sim_outputs = {"Date": [], "Mean": []}
    df = pd.read_csv(sim_output_filepath + "/new_infections_mean.csv")
    current_date = start_date
    for week in range(0, len(df)):
        incidence = (df.iloc[week, region + 1]) / pop_size * 100000
        sim_outputs["Date"].append(current_date)
        sim_outputs["Mean"].append(incidence)
        current_date += pd.Timedelta(days=7)
    sim_outputs = pd.DataFrame(sim_outputs)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            marker='o', linestyle='-', color=colors['dark red'], label='simulated')

    if (real_data_file != ""):
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        filtered_df = real_df[(real_df["Datum"] >= start_date)
                              & (real_df["Datum"] <= end_date) & (real_df["Region"] == region_name)]
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
                f"incidence_{region_name}.png")


def plot_infected_mean_std(sim_output_filepath, save_dir, start_date, region_name, tmax, plot_new_inf, pop_size = 100000, col="C2", col_var="M020", real_data_file = "", real_Altersgruppe = ""):
    sim_outputs = {"Date": [], "Mean": [], "Std": []}
    df_mean = pd.read_csv(sim_output_filepath + "/means.csv")
    df_moments = pd.read_csv(sim_output_filepath + "/moments.csv")
    current_date = start_date

    for t in range(0, len(df_mean)):
        sim_outputs["Date"].append(current_date + pd.to_timedelta(t/10., unit="D"))
        sim_outputs["Mean"].append(df_mean[col].iloc[t])
        sim_outputs["Std"].append(np.sqrt(df_moments[col_var].iloc[t]))

    sim_outputs = pd.DataFrame(sim_outputs)
    interval = int(8*7/0.1)

    # Plot mean +- std
    fig, ax = plt.subplots(figsize=(12, 6))
        
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            color=colors['middle blue'], label=r'$\mu_I$', zorder=20)
    ax.fill_between(sim_outputs["Date"], sim_outputs["Mean"] - 2*sim_outputs["Std"], sim_outputs["Mean"] + 2*sim_outputs["Std"],
                    color=colors['middle blue'], alpha=0.5, label=r'$\mu_I \pm 2\sigma_I$', zorder=5)
    ax.set_xlabel("Date")
    ax.set_ylabel("Infected [#]")
    ax.set_xticks(sim_outputs["Date"].iloc[::interval])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    
    if (real_data_file != ""):
        ax1 = ax.twinx()
        if(plot_new_inf):
            sim_outputs_new_inf = {"Date": [], "Mean": []}
            df = pd.read_csv(sim_output_filepath + "/new_infections_mean.csv")
            current_date = start_date
            for week in range(0, len(df)):
                incidence = (df.iloc[week, region + 1]) / pop_size * 100000
                sim_outputs_new_inf["Date"].append(current_date)
                sim_outputs_new_inf["Mean"].append(incidence)
                current_date += pd.Timedelta(days=7)
            sim_outputs_new_inf = pd.DataFrame(sim_outputs_new_inf)
            ax1.plot(sim_outputs_new_inf["Date"], sim_outputs_new_inf["Mean"],
            marker='^', linestyle='-', color=colors['dark red'], label='simulated incidence', zorder=10)
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        filtered_df = real_df[(real_df["Datum"] >= start_date)
                              & (real_df["Datum"] <= end_date) & (real_df["Region"] == region_name) & (real_df["Altersgruppe"] == real_Altersgruppe)]
        ax1.scatter(filtered_df["Datum"], filtered_df["Inzidenz"],
                   marker='x', color='black', label='real incidence', zorder=10)
        ax1.set_ylabel("Real incidence")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax1.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    else:
        lines, labels = ax.get_legend_handles_labels()
        ax.legend(lines, labels, loc='upper left')
        
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir + f"mean_std_{region_name}.png")

def plot_infected_percentiles(sim_output_filepath, save_dir, start_date, p_lower, p_upper, region_name, tmax, col="C2", real_data_file = "", real_Altersgruppe = ""):
    sim_outputs = {"Date": [], "Mean": [], "Lower": [], "Upper": []}
    df_mean = pd.read_csv(sim_output_filepath + "/means.csv")
    df_lower = pd.read_csv(sim_output_filepath + f"/{p_lower}.csv")
    df_upper = pd.read_csv(sim_output_filepath + f"/{p_upper}.csv")
    current_date = start_date

    for t in range(0, len(df_mean)):
        sim_outputs["Date"].append(current_date + pd.to_timedelta(t/10., unit="D"))
        sim_outputs["Mean"].append(df_mean[col].iloc[t])
        sim_outputs["Lower"].append(df_lower[col].iloc[t])
        sim_outputs["Upper"].append(df_upper[col].iloc[t])

    sim_outputs = pd.DataFrame(sim_outputs)
    interval = int(8*7/0.1)

    # Plot mean and percentiles
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(sim_outputs["Date"], sim_outputs["Mean"],
            color=colors['middle blue'], label=r'$\mu_I$')
    ax.fill_between(sim_outputs["Date"], sim_outputs["Lower"], sim_outputs["Upper"],
                    color=colors['middle blue'], alpha=0.5, label=r'$90$-quantile')
    ax.set_xlabel("Date")
    ax.set_ylabel("Infected [#]")
    ax.set_xticks(sim_outputs["Date"].iloc[::interval])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    
    if (real_data_file != ""):
        ax1 = ax.twinx()
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        filtered_df = real_df[(real_df["Datum"] >= start_date)
                              & (real_df["Datum"] <= end_date) & (real_df["Region"] == region_name) & (real_df["Altersgruppe"] == real_Altersgruppe)]
        ax1.scatter(filtered_df["Datum"], filtered_df["Inzidenz"],
                   marker='x', color='black', label='real')
        ax1.set_ylabel("Real incidence")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax1.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    else:
        lines, labels = ax.get_legend_handles_labels()
        ax.legend(lines, labels, loc='upper left')
        
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir + f"percentiles_{region_name}.png")

sim_output_dir = "V:/bick_ju/TemporalHybrid/SMM/" #"V:/bick_ju/TemporalHybrid/Moments/config_influenza_germany/truncation/closure_order_3/0.000000/"#"V:/bick_ju/TemporalHybrid/SMM/"
grippe_web_data = "C:/Users/bick_ju/Documents/repos/grippeweb_data/ILI_df_Germany.csv"
start_date = pd.to_datetime("2016-08-01")
save_dir = "H:/Documents/TemporalHybridModel/SMM/"#"H:/Documents/TemporalHybridModel/Moments/"#"H:/Documents/TemporalHybridModel/SMM/"
config = "config_influenza_germany"
sim_output_dir += config
save_dir += config + "/"
region = 0
pop_size = 100000
region_name = "Bundesweit"
tmax = 3 * 365 - 3

os.makedirs(save_dir, exist_ok=True)


plot_infected_mean_std(sim_output_dir,
              save_dir, start_date, region_name, tmax, plot_new_inf=True, real_data_file=grippe_web_data, real_Altersgruppe="00+")

plot_infected_percentiles(sim_output_dir,
              save_dir, start_date, "p05", "p95", region_name, tmax, real_data_file=grippe_web_data, real_Altersgruppe="00+")

plot_new_infections(sim_output_dir, save_dir,
                    grippe_web_data, start_date, tmax, region, pop_size, region_name)
