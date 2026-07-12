from _simulation_stochastic import run_germany
import pandas as pd

# 0 is setup with pop size 100000 and 1 with posize corresponding to full Germany
full_or_scaled = 0
num_runs = 1

trams_rate = 0.0000018
gamma = 1./7.
nu = 1./14.
I0 = 1000
R0 = 20000
peaks = [-10, -5, -5]
rhos = [0.5, 0.7, 0.7]
sigmas = [100, 100, 100]

res = run_germany(full_or_scaled, num_runs, trams_rate,
                  gamma, nu, I0, R0, peaks, rhos, sigmas)

target_file = "/Users/julia/repos/grippeweb_data/ILI_df_Germany.csv"

start_date = pd.to_datetime("2016-08-01")
num_days = 3 * 365 - 4

real_df = pd.read_csv(target_file, parse_dates=["Datum"])
end_date = start_date + pd.Timedelta(days=num_days)
filtered_df = real_df[(real_df["Datum"] >= start_date)
                      & (real_df["Datum"] <= end_date)]

target = list(filtered_df.Inzidenz)

print()
