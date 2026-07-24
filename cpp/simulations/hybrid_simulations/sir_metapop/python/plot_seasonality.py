from settings import *

beta0s = [1, 1]
rhos = [[0.5, 0.5, 0.5], [0.3, 0.5, 0.7]]
peaks = [[182, 182, 182], [160, 210, 150]]
sigmas = [[50, 50, 50], [20, 50, 100]]

def gaussian(t, t_peak_season, sigma):    
    return 1. / (sigma * sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.pow((t - t_peak_season) / sigma, 2))

def zeta(t, season, start_day):
    tau_minus = start_day + season * 365.
    tau_plus = tau_minus + 365.
    if (t <= tau_minus + 30):
        return (t - tau_minus + 30.) / 60.
    elif(t >= tau_plus - 30. and t <= tau_plus):
        return (-t + tau_plus + 30.) / 60.
    return 1

def calculate_seasonality_factor(t, peaks, rhos, sigmas, start_day):
    season           = int((t - start_day) /
                         365.)
    t_peak_season = season * 365. + peaks[season]
    rho = rhos[season]
    sigma = sigmas[season]
    return rho + (1 - rho) * gaussian(t, t_peak_season, sigma) / gaussian(0, 0, sigma) * zeta(t, season, start_day)

fig, ax = plt.subplots(figsize=(5, 3))
for setting in range(len(rhos)):
    for season in range(len(rhos[setting])):
        first_season_start_day = 365 * season
        time = [t for t in range(first_season_start_day + 365)]
        y = []
        for t in time:
            y.append(calculate_seasonality_factor(t-365, peaks[setting], rhos[setting], sigmas[setting], first_season_start_day))
        ax.plot(time, y)
ax.set_xlabel("Time [days]")
ax.set_ylabel("Transmission rate")
fig.savefig("Seasonality.png")