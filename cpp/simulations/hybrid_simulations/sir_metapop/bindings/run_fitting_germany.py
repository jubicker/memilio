from _simulation_stochastic import run_germany
import pandas as pd
import pyabc
import numpy as np

# Load fitting data
target_file = "ILI_df_Germany.csv"

start_date = pd.to_datetime("2016-08-01")
num_days = 3 * 365 - 4

real_df = pd.read_csv(target_file, parse_dates=["Datum"])
end_date = start_date + pd.Timedelta(days=num_days)
filtered_df = real_df[(real_df["Datum"] >= start_date)
                      & (real_df["Datum"] <= end_date)]
target = np.array(filtered_df.Inzidenz)
obs_data = dict(data=target)

#setup parameters
full_or_scaled = 0
num_runs = 1

#fixed parameters that are not fitted:
gamma = 1./7.
nu = 1./14.

# Model function that gets parameters, runs simulation and returns result dictionary
def model(parameters):
    res = run_germany(
        full_or_scaled, 
        num_runs, 
        10**(-1*parameters["trams_rate"]),
        gamma, 
        nu, 
        10**parameters["I0"], 
        10**parameters["R0"], 
        [round(parameters["peaks_year_1"]),
        round(parameters["peaks_year_2"]),
        round(parameters["peaks_year_3"])], 
        [parameters["rhos_year_1"],
        parameters["rhos_year_2"],
        parameters["rhos_year_3"]], 
        [parameters["sigmas_year_1"],
        parameters["sigmas_year_2"],
        parameters["sigmas_year_3"]])
    # We only use one region, so we can directly simplify the data
    return{"data": np.array([week[0] for week in res])}

# Define the prior for all other parameters
prior = pyabc.Distribution(
    trams_rate=pyabc.RV("uniform", 3, 7),
    I0=pyabc.RV("uniform", 2, 4),
    R0=pyabc.RV("uniform", 2, 4),
    peaks_year_1=pyabc.RV("norm", loc = 0, scale = 20),
    peaks_year_2=pyabc.RV("norm", loc = 0, scale = 20),
    peaks_year_3=pyabc.RV("norm", loc = 0, scale = 20), 
    rhos_year_1=pyabc.RV("uniform", 0, 1),
    rhos_year_2=pyabc.RV("uniform", 0, 1),
    rhos_year_3=pyabc.RV("uniform", 0, 1), 
    sigmas_year_1=pyabc.RV("uniform", 0, 200),
    sigmas_year_2=pyabc.RV("uniform", 0, 200),
    sigmas_year_3=pyabc.RV("uniform", 0, 200),
)

# Define the fitting problem
abc = pyabc.ABCSMC(model, prior, pyabc.distance.PNormDistance(), population_size = 100)

# Create a database
db_path = "sqlite:///influenca_fitting.db"
abc.new(db_path, obs_data)

# Run the fitting
history = abc.run(max_nr_populations=10, minimum_epsilon=0.1, min_acceptance_rate=0.01)





# trams_rate = 0.0000018
# I0 = 1000
# R0 = 20000
# peaks = [-10, -5, -5]
# rhos = [0.5, 0.7, 0.7]
# sigmas = [100, 100, 100]
# print()
