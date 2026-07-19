from _simulation_stochastic import run_germany
import pandas as pd
import pyabc
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import os

# Load fitting data
target_file = "/Users/julia/repos/grippeweb_data/ILI_df_Germany.csv"

start_date = pd.to_datetime("2016-08-01")
num_days = 3 * 365 - 4

real_df = pd.read_csv(target_file, parse_dates=["Datum"])
end_date = start_date + pd.Timedelta(days=num_days)
filtered_df = real_df[(real_df["Datum"] >= start_date)
                      & (real_df["Datum"] <= end_date)]
target = np.array(filtered_df.Inzidenz)
obs_data = dict(data=target)

# setup parameters
full_or_scaled = 0
use_new_inf = True
num_runs = 1

# fixed parameters that are not fitted:
gamma = 1./7.

# Model function that gets parameters, runs simulation and returns result dictionary


def model(parameters):
    res = run_germany(
        full_or_scaled, use_new_inf,
        num_runs,
        10**(-1*parameters["trams_rate"]),
        gamma,
        1./parameters["T_R"],
        parameters["I0"],
        parameters["R0"],
        [parameters["peaks_year_1"],
         parameters["peaks_year_2"],
         parameters["peaks_year_3"]],
        [parameters["rhos_year_1"],
         parameters["rhos_year_2"],
         parameters["rhos_year_3"]],
        [parameters["sigmas_year_1"],
         parameters["sigmas_year_2"],
         parameters["sigmas_year_3"]])
    # We only use one region, so we can directly simplify the data
    return {"data": np.array([week[0] for week in res])}


domain_I0 = np.arange(0, 2000)
domain_R0 = np.arange(0, 30000)
domain_peaks = np.arange(-30, 30)

# Define the prior for all other parameters
# prior = pyabc.Distribution(
#     trams_rate=pyabc.RV("uniform", 4, 2),
#     I0=pyabc.RV("rv_discrete", values=(domain_I0, [1 / 2000] * 2000)),
#     R0=pyabc.RV("rv_discrete", values=(domain_R0, [1 / 30000] * 30000)),
#     peaks_year_1=pyabc.RV("rv_discrete", values=(domain_peaks, [1 / 60] * 60)),
#     peaks_year_2=pyabc.RV("rv_discrete", values=(domain_peaks, [1 / 60] * 60)),
#     peaks_year_3=pyabc.RV("rv_discrete", values=(domain_peaks, [1 / 60] * 60)),
#     rhos_year_1=pyabc.RV("uniform", 0, 1),
#     rhos_year_2=pyabc.RV("uniform", 0, 1),
#     rhos_year_3=pyabc.RV("uniform", 0, 1),
#     sigmas_year_1=pyabc.RV("uniform", 0, 200),
#     sigmas_year_2=pyabc.RV("uniform", 0, 200),
#     sigmas_year_3=pyabc.RV("uniform", 0, 200),
# )

prior = pyabc.Distribution(
    trams_rate=pyabc.RV("uniform", 4, 2),
    T_R=pyabc.RV("uniform", 10, 355),
    I0=pyabc.RV("uniform", 0, 3000),
    R0=pyabc.RV("uniform", 0, 30000),
    peaks_year_1=pyabc.RV("norm", loc=0, scale=20),
    peaks_year_2=pyabc.RV("norm", loc=0, scale=20),
    peaks_year_3=pyabc.RV("norm", loc=0, scale=20),
    rhos_year_1=pyabc.RV("uniform", 0, 1),
    rhos_year_2=pyabc.RV("uniform", 0, 1),
    rhos_year_3=pyabc.RV("uniform", 0, 1),
    sigmas_year_1=pyabc.RV("uniform", 0, 200),
    sigmas_year_2=pyabc.RV("uniform", 0, 200),
    sigmas_year_3=pyabc.RV("uniform", 0, 200),
)


transition = pyabc.AggregatedTransition(
    mapping={
        'trams_rate': pyabc.MultivariateNormalTransition(),
        'T_R': pyabc.MultivariateNormalTransition(),
        'I0': pyabc.DiscreteJumpTransition(domain=domain_I0, p_stay=0.7),
        'R0': pyabc.DiscreteJumpTransition(domain=domain_R0, p_stay=0.7),
        'peaks_year_1': pyabc.DiscreteJumpTransition(domain=domain_peaks, p_stay=0.7),
        'peaks_year_2': pyabc.DiscreteJumpTransition(domain=domain_peaks, p_stay=0.7),
        'peaks_year_3': pyabc.DiscreteJumpTransition(domain=domain_peaks, p_stay=0.7),
        'rhos_year_1': pyabc.MultivariateNormalTransition(),
        'rhos_year_2': pyabc.MultivariateNormalTransition(),
        'rhos_year_3': pyabc.MultivariateNormalTransition(),
        'sigmas_year_1': pyabc.MultivariateNormalTransition(),
        'sigmas_year_2': pyabc.MultivariateNormalTransition(),
        'sigmas_year_3': pyabc.MultivariateNormalTransition()
    }
)


def distance(x, x0):
    return np.square(x["data"] - x0["data"]).mean()


def weighted_quantiles(simulations, weights, qs):
    """
    Compute weighted quantiles column-wise.

    Parameters
    ----------
    simulations : array-like, shape (n, m)
        n simulations, m time points (or dimensions).
    weights : array-like, shape (n,)
        Nonnegative weights.
    qs : float or array-like
        Quantile(s) in [0, 1]. Example: 0.5 or (0.25, 0.5, 0.75).

    Returns
    -------
    quantiles : ndarray
        If qs is scalar -> shape (m,)
        If qs is iterable of length k -> shape (k, m)
    """
    x = np.asarray(simulations)
    w = np.asarray(weights)
    if x.ndim != 2:
        raise ValueError("simulations must be 2D (n, m)")
    if w.ndim != 1 or w.shape[0] != x.shape[0]:
        raise ValueError("weights must be 1D with length n")
    if np.any(w < 0):
        raise ValueError("weights must be nonnegative")
    tot = w.sum()
    if tot <= 0:
        raise ValueError("sum of weights must be > 0")
    qs = np.atleast_1d(qs)
    # sort each column
    idx = np.argsort(x, axis=0)                   # (n, m)
    x_sorted = np.take_along_axis(x, idx, axis=0)
    # reorder weights the same way (no huge broadcasted temp arrays)
    w_sorted = w[idx]                             # (n, m)
    cw = np.cumsum(w_sorted, axis=0)              # (n, m)
    results = []
    for q in qs:
        cutoff = q * tot
        k = (cw >= cutoff).argmax(axis=0)         # (m,)
        results.append(x_sorted[k, np.arange(x.shape[1])])
    result = np.vstack(results)
    if result.shape[0] == 1:
        return result[0]                          # (m,)
    return result


def plot_new_infections_cis(sim_output_matrix, save_dir, real_data_file, start_date, tmax, region, pop_size, region_name):
    # Median
    median_curve = {"Date": [], "Mean": []}
    current_date = start_date
    for week in range(0, len(sim_output_matrix[1])):
        incidence = (sim_output_matrix[1][week]) / pop_size * 100000
        median_curve["Date"].append(current_date)
        median_curve["Mean"].append(incidence)
        current_date += pd.Timedelta(days=7)
    median_curve = pd.DataFrame(median_curve)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(median_curve["Date"], median_curve["Mean"],
            marker='o', linestyle='-', color='cornflowerblue', label='simulated median')
    # Lower Bound
    lower_curve = {"Date": [], "Mean": []}
    current_date = start_date
    for week in range(0, len(sim_output_matrix[0])):
        incidence = (sim_output_matrix[0][week]) / pop_size * 100000
        lower_curve["Date"].append(current_date)
        lower_curve["Mean"].append(incidence)
        current_date += pd.Timedelta(days=7)
    lower_curve = pd.DataFrame(lower_curve)
    # Upper Bound
    upper_curve = {"Date": [], "Mean": []}
    current_date = start_date
    for week in range(0, len(sim_output_matrix[2])):
        incidence = (sim_output_matrix[2][week]) / pop_size * 100000
        upper_curve["Date"].append(current_date)
        upper_curve["Mean"].append(incidence)
        current_date += pd.Timedelta(days=7)
    upper_curve = pd.DataFrame(upper_curve)
    ax.fill_between(upper_curve["Date"], lower_curve["Mean"], upper_curve["Mean"],
                    color='cornflowerblue', alpha=0.3, label='CrI')
    if (real_data_file != ""):
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        filtered_df = real_df[(real_df["Datum"] >= start_date)
                              & (real_df["Datum"] <= end_date) & (real_df["Region"] == region_name)]
        ax.scatter(filtered_df["Datum"], filtered_df["Inzidenz"],
                   marker='x', color='black', label='real')
    ax.set_xlabel("Date")
    ax.set_ylabel("Incidence")
    ax.set_xticks(median_curve["Date"].iloc[::8])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, loc='upper left')
    plt.grid()
    plt.tight_layout()
    fig.savefig(save_dir +
                f"incidence_{region_name}.png")


if __name__ == "__main__":
    dir_path = "/Users/julia/repos/fork/memilio/cpp/simulations/hybrid_simulations/sir_metapop/bindings/output/"
    os.makedirs(dir_path, exist_ok=True)
    # Create a database
    db_path = "sqlite:///" + dir_path + "influenca_fitting3.db"
    run = True
    if (run):
        # Define the fitting problem
        abc = pyabc.ABCSMC(model, prior, distance, population_size=pyabc.populationstrategy.AdaptivePopulationSize(
            300, mean_cv=0.05, max_population_size=1000))

        abc.new(db_path, obs_data)

        # Run the fitting
        history = abc.run(max_nr_populations=10, minimum_epsilon=0.1)

    else:
        history = pyabc.History(db_path, create=False)

    pop = history.get_population()
    best_particle = min(pop.particles, key=lambda p: p.distance)
    output = pd.DataFrame({"Mean": best_particle.sum_stat["data"]})
    output.to_csv(dir_path + "new_infections_mean.csv")

    top10_particles = sorted(pop.particles, key=lambda p: p.distance)[:10]
    top10_df = pd.DataFrame([dict(p.parameter) for p in top10_particles])
    top10_df.insert(0, "distance", [p.distance for p in top10_particles])
    top10_df.to_csv(dir_path + "top10_particles.csv", index=False)

    df, w = history.get_distribution()

    pyabc.visualization.plot_kde_matrix(df, w)
    plt.savefig(dir_path + "posterior_kde_matrix.png")

    _fig, _arr_ax = plt.subplots(1, 5, figsize=(10, 2.5))
    _arr_ax = _arr_ax.flatten()
    pyabc.visualization.plot_sample_numbers(history, ax=_arr_ax[0])
    _arr_ax[0].get_legend().remove()
    pyabc.visualization.plot_walltime(history, ax=_arr_ax[1], unit='h')
    _arr_ax[1].get_legend().remove()
    pyabc.visualization.plot_epsilons(history, ax=_arr_ax[2])
    pyabc.visualization.plot_effective_sample_sizes(history, ax=_arr_ax[3])
    pyabc.visualization.plot_acceptance_rates_trajectory(
        history, ax=_arr_ax[4])
    plt.savefig(dir_path + "stats.png")

    result_matrix = np.vstack([x.sum_stat["data"] for x in pop.particles])
    weighted_results = weighted_quantiles(result_matrix, w, (0.05, 0.5, 0.95))

    plot_new_infections_cis(weighted_results, dir_path, target_file, pd.Timestamp(
        "2016-08-01"), 3*365-4, 0, 100000, "Bundesweit")

    # trams_rate = 0.0000018
    # I0 = 1000
    # R0 = 20000
    # peaks = [-10, -5, -5]
    # rhos = [0.5, 0.7, 0.7]
    # sigmas = [100, 100, 100]
    # print()
