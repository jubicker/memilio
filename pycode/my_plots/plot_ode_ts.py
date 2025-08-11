from settings import *

result_file = sim_results + 'ODE/sir.csv'

df = pd.read_csv(result_file)
comps = list(df.columns)
comps.remove('Time')

width_in_cm = 7
height_in_cm = 5.5

figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)
for c in comps:
    ax.plot(df.Time, df[c], color=colors[c], label=c)
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder+'ode_sir.png', dpi=dpi)

