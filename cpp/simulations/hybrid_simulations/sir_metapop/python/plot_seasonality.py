from settings import *

beta0s = [1, 1]
rhos = [[0.5, 0.5, 0.5], [0.3, 0.5, 0.7]]
peaks = [[182, 182, 182], [160, 210, 150]]
sigmas = [[40, 40, 40], [15, 40, 60]]

setting_colors = {0: colors["brown"], 1: colors["teal"]}

def gaussian(t, t_peak_season, sigma):    
    return 1. / (sigma * sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.pow((t - t_peak_season) / sigma, 2))

def zeta(t, start_day):
    tau_minus = start_day
    tau_plus = tau_minus + 365.
    if (t <= tau_minus + 30):
        return (t - tau_minus + 30.) / 60.
    elif(t >= tau_plus - 30. and t <= tau_plus):
        return (-t + tau_plus + 30.) / 60.
    return 1

def calculate_seasonality_factor(t, peak, rho, sigma):
    return rho + (1 - rho) * gaussian(t, peak, sigma) / gaussian(0, 0, sigma) * zeta(t, 0)

fig, ax = plt.subplots(figsize=(2.8, 1.6), dpi=dpi)
for setting in range(len(rhos)):
    for season in range(len(rhos[setting])):
        first_season_start_day = 365 * season
        time = [t for t in range(first_season_start_day, first_season_start_day + 365)]
        y = []
        for t in time:
            y.append(calculate_seasonality_factor(t-first_season_start_day, peaks[setting][season], rhos[setting][season], sigmas[setting][season]))
        ax.plot(time, y, color=setting_colors[setting])
ax.set_xlabel("Time [days]")
ax.set_ylabel(r"$\lambda$ [normalized]")
fig.subplots_adjust(left=0.2, bottom=0.27, top=0.98, right=0.98)
fig.savefig("Seasonality.png", dpi=dpi)