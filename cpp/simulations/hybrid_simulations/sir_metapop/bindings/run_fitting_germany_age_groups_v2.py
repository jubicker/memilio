from _simulation_stochastic import run_age_groups
import pandas as pd
import pyabc
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import os

# Load fitting data
target_file = "/p/project1/loki/bicker1/memilio/ILI_df_age_groups.csv"
age_group_to_index = {"0-4": 0, "5-14": 1, "15-34": 2, "35-59": 3, "60+": 4}

start_date = pd.to_datetime("2016-08-01")
num_days = 3 * 365 - 4

real_df = pd.read_csv(target_file, parse_dates=["Datum"])
end_date = start_date + pd.Timedelta(days=num_days)
target = []  # For every time point one list with one value for each age group
filtered_df = real_df[(real_df["Datum"] >= start_date)
                      & (real_df["Datum"] <= end_date)]
for date in filtered_df.Datum.unique():
    ag_values = []
    for ag in age_group_to_index.keys():
        value = filtered_df[(filtered_df.Datum == date) &
                            (filtered_df.Altersgruppe == ag)]
        ag_values.append(value.Inzidenz.iloc[0])
    target.append(ag_values)

target_arr = np.array(target)          # shape (T, 5)
group_scale = target_arr.std(axis=0)   # per-age-group std, shape (5,)
obs_data = dict(data=target)

# setup parameters
full_or_scaled = 0
use_new_inf = True
num_runs = 1

# fixed parameters that are not fitted:
gamma = 1./7.
nu = 1./149.
peaks = [21, 37, 24]
sigmas = [10, 15, 5]
rhos = [0.79, 0.78, 0.77]

# Model function that gets parameters, runs simulation and returns result dictionary


def model(parameters):
    res = run_age_groups(
        full_or_scaled, use_new_inf,
        num_runs,
        [10**(-1*parameters["trams_rate_group_1"]), 10**(-1*parameters["trams_rate_group_2"]), 10**(-1 *
                                                                                                    parameters["trams_rate_group_3"]), 10**(-1*parameters["trams_rate_group_4"]), 10**(-1*parameters["trams_rate_group_5"])],
        [[parameters["f11"], parameters["f12"], parameters["f13"], parameters["f14"], parameters["f15"]],
         [parameters["f21"], parameters["f22"], parameters["f23"], parameters["f24"], parameters["f25"]],
         [parameters["f31"], parameters["f32"], parameters["f33"], parameters["f34"], parameters["f35"]],
         [parameters["f41"], parameters["f42"], parameters["f43"], parameters["f44"], parameters["f45"]],
         [parameters["f51"], parameters["f52"], parameters["f53"], parameters["f54"], parameters["f55"]]],
        gamma,
        nu,
        [parameters["I0_group_1"], parameters["I0_group_2"], parameters["I0_group_3"],
            parameters["I0_group_4"], parameters["I0_group_5"]],
        [parameters["R0_group_1"], parameters["R0_group_2"], parameters["R0_group_3"],
            parameters["R0_group_4"], parameters["R0_group_5"]],
        peaks,
        rhos,
        sigmas)
    return {"data": res}


prior = pyabc.Distribution(
    trams_rate_group_1=pyabc.RV("uniform", 4, 3),
    trams_rate_group_2=pyabc.RV("uniform", 4, 3),
    trams_rate_group_3=pyabc.RV("uniform", 4, 3),
    trams_rate_group_4=pyabc.RV("uniform", 4, 3),
    trams_rate_group_5=pyabc.RV("uniform", 4, 3),
    f11 = pyabc.RV("uniform", 0, 1),
    f12 = pyabc.RV("uniform", 0, 1),
    f13 = pyabc.RV("uniform", 0, 1),
    f14 = pyabc.RV("uniform", 0, 1),
    f15 = pyabc.RV("uniform", 0, 1),
    f21 = pyabc.RV("uniform", 0, 1),
    f22 = pyabc.RV("uniform", 0, 1),
    f23 = pyabc.RV("uniform", 0, 1),
    f24 = pyabc.RV("uniform", 0, 1),
    f25 = pyabc.RV("uniform", 0, 1),
    f31 = pyabc.RV("uniform", 0, 1),
    f32 = pyabc.RV("uniform", 0, 1),
    f33 = pyabc.RV("uniform", 0, 1),
    f34 = pyabc.RV("uniform", 0, 1),
    f35 = pyabc.RV("uniform", 0, 1),
    f41 = pyabc.RV("uniform", 0, 1),
    f42 = pyabc.RV("uniform", 0, 1),
    f43 = pyabc.RV("uniform", 0, 1),
    f44 = pyabc.RV("uniform", 0, 1),
    f45 = pyabc.RV("uniform", 0, 1),
    f51 = pyabc.RV("uniform", 0, 1),
    f52 = pyabc.RV("uniform", 0, 1),
    f53 = pyabc.RV("uniform", 0, 1),
    f54 = pyabc.RV("uniform", 0, 1),
    f55 = pyabc.RV("uniform", 0, 1),
    I0_group_1=pyabc.RV("uniform", 0, 3000),
    I0_group_2=pyabc.RV("uniform", 0, 3000),
    I0_group_3=pyabc.RV("uniform", 0, 3000),
    I0_group_4=pyabc.RV("uniform", 0, 3000),
    I0_group_5=pyabc.RV("uniform", 0, 3000),
    R0_group_1=pyabc.RV("uniform", 0, 30000),
    R0_group_2=pyabc.RV("uniform", 0, 30000),
    R0_group_3=pyabc.RV("uniform", 0, 30000),
    R0_group_4=pyabc.RV("uniform", 0, 30000),
    R0_group_5=pyabc.RV("uniform", 0, 30000)
)


def distance(x, x0):
    sim = np.asarray(x["data"])
    obs = np.asarray(x0["data"])
    return np.mean(((sim - obs)) ** 2)


def weighted_quantiles(simulations, weights, qs):
    """
    Compute weighted quantiles across simulations (first axis), independently
    for every entry of the remaining dimensions.

    Parameters
    ----------
    simulations : array-like, shape (n, ...)
        n simulations, followed by arbitrary output dimensions, e.g.
        (n, m) for m time points, or (n, m, a) for m time points and
        a age groups.
    weights : array-like, shape (n,)
        Nonnegative weights.
    qs : float or array-like
        Quantile(s) in [0, 1]. Example: 0.5 or (0.25, 0.5, 0.75).

    Returns
    -------
    quantiles : ndarray
        If qs is scalar -> shape simulations.shape[1:]
        If qs is iterable of length k -> shape (k,) + simulations.shape[1:]
    """
    x = np.asarray(simulations)
    w = np.asarray(weights)
    if x.ndim < 2:
        raise ValueError("simulations must have at least 2 dimensions (n, ...)")
    if w.ndim != 1 or w.shape[0] != x.shape[0]:
        raise ValueError("weights must be 1D with length n")
    if np.any(w < 0):
        raise ValueError("weights must be nonnegative")
    tot = w.sum()
    if tot <= 0:
        raise ValueError("sum of weights must be > 0")
    qs = np.atleast_1d(qs)

    n = x.shape[0]
    out_shape = x.shape[1:]
    x2 = x.reshape(n, -1)                          # (n, m) with m = prod(out_shape)

    # sort each column
    idx = np.argsort(x2, axis=0)                   # (n, m)
    x_sorted = np.take_along_axis(x2, idx, axis=0)
    # reorder weights the same way (no huge broadcasted temp arrays)
    w_sorted = w[idx]                              # (n, m)
    cw = np.cumsum(w_sorted, axis=0)                # (n, m)
    results = []
    for q in qs:
        cutoff = q * tot
        k = (cw >= cutoff).argmax(axis=0)          # (m,)
        results.append(x_sorted[k, np.arange(x2.shape[1])].reshape(out_shape))
    result = np.stack(results, axis=0)
    if result.shape[0] == 1:
        return result[0]                           # out_shape
    return result                                  # (k,) + out_shape


def plot_new_infections_cis(sim_output_matrix, save_dir, real_data_file, start_date, tmax, region, pop_size, region_name, age_group_to_index):
    """
    Plot simulated median + credible interval (and real data, if given) of
    new infections per 100,000 population, with one subplot per age group.

    Parameters
    ----------
    sim_output_matrix : ndarray, shape (3, num_weeks, num_age_groups)
        Lower bound (index 0), median (index 1) and upper bound (index 2) of
        the credible interval, as returned by weighted_quantiles with
        qs=(low, 0.5, high).
    pop_size : scalar or array-like of length num_age_groups
        Population size used to normalize incidence. A scalar is broadcast
        to all age groups.
    age_group_to_index : dict
        Maps age group label (e.g. "0-4") to its column index in
        sim_output_matrix, in the same order used when assembling the
        fitting target.
    """
    num_age_groups = sim_output_matrix.shape[2]
    pop_size = np.broadcast_to(
        np.asarray(pop_size, dtype=float), (num_age_groups,))

    dates = [start_date + pd.Timedelta(days=7 * week)
             for week in range(sim_output_matrix.shape[1])]

    real_df = None
    if real_data_file != "":
        real_df = pd.read_csv(real_data_file, parse_dates=["Datum"])
        end_date = start_date + pd.Timedelta(days=tmax)
        real_df = real_df[(real_df["Datum"] >= start_date)
                          & (real_df["Datum"] <= end_date) & (real_df["Region"] == region_name)]

    ncols = 3
    nrows = int(np.ceil(num_age_groups / ncols))
    fig, _ = plt.subplots(
        nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False)
    axes = fig.axes

    for ag, idx in age_group_to_index.items():
        ax = axes[idx]
        lower = sim_output_matrix[0, :, idx] / pop_size[idx] * 100000
        median = sim_output_matrix[1, :, idx] / pop_size[idx] * 100000
        upper = sim_output_matrix[2, :, idx] / pop_size[idx] * 100000

        ax.plot(dates, median, marker='o', linestyle='-',
                color='cornflowerblue', label='simulated median')
        ax.fill_between(dates, lower, upper, color='cornflowerblue',
                        alpha=0.3, label='CrI')

        if real_df is not None:
            ag_df = real_df[real_df["Altersgruppe"] == ag]
            ax.scatter(ag_df["Datum"], ag_df["Inzidenz"],
                       marker='x', color='black', label='real')

        ax.set_title(f"Age group {ag}")
        ax.set_xlabel("Date")
        ax.set_ylabel("Incidence")
        ax.set_xticks(dates[::8])
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        ax.grid()

    for idx in range(num_age_groups, len(axes)):
        axes[idx].set_visible(False)

    lines, labels = axes[0].get_legend_handles_labels()
    fig.legend(lines, labels, loc='upper center', ncol=3)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(save_dir + f"incidence_{region_name}_age_groups.png")


if __name__ == "__main__":
    dir_path = "/p/project1/loki/bicker1/memilio/cpp/simulations/hybrid_simulations/sir_metapop/bindings/output_ag/"
    os.makedirs(dir_path, exist_ok=True)
    # Create a database
    db_path = "sqlite:///" + dir_path + "influenca_fitting3.db"
    run = True
    load = False
    if (run):
        # Define the fitting problem
        abc = pyabc.ABCSMC(model, prior, distance, population_size=4000) #pyabc.populationstrategy.AdaptivePopulationSize( 500, mean_cv=0.05, max_population_size=3000)
        
        if(load):
            abc.load(db_path, 1)
        else:
            abc.new(db_path, obs_data)

        # Run the fitting
        history = abc.run(max_nr_populations=50, minimum_epsilon=0.1, min_acceptance_rate=0.01)

    else:
        history = pyabc.History(db_path, create=False)

    pop = history.get_population()
    best_particle = min(pop.particles, key=lambda p: p.distance)
    output = pd.DataFrame({"Mean": [best_particle.sum_stat["data"][t, :] for t in range(best_particle.sum_stat["data"].shape[0])]})
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

    result_matrix = np.stack(
        [np.asarray(x.sum_stat["data"]) for x in pop.particles])  # (n, T, A)
    weighted_results = weighted_quantiles(
        result_matrix, w, (0.05, 0.5, 0.95))  # (3, T, A)

    plot_new_infections_cis(weighted_results, dir_path, target_file, pd.Timestamp(
        "2016-08-01"), 3*365-4, 0, 100000, "Bundesweit", age_group_to_index)
