from settings import *

num_sims = 500
dfs = []
for i in range(num_sims):
    result_file = sim_results + f'ABM/{i}_comps.csv'
    df= pd.read_csv(result_file)
    df = df.drop(columns=['I_NS', 'I_Sy', 'I_Sev', 'I_Crit', 'D'])
    df = df.rename(columns={'E': 'I'})
    dfs.append(df)
    
# Assumes all have the same time points in the same order
time = dfs[0]['Time']
comps = list(dfs[0].columns)
comps.remove('Time')

width_in_cm = 7
height_in_cm = 5.5

figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)
for c in comps:
    matrix = np.array([df[c].values for df in dfs])
    mean_val = np.mean(matrix, axis=0)
    p5 = np.percentile(matrix, 5, axis=0)
    p95 = np.percentile(matrix, 95, axis=0)
    ax.plot(time, mean_val, color=colors[c], label=c)
    ax.fill_between(time, p5, p95, color=colors[c], alpha=0.2, lw=0)
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + 'abm_sir.png', dpi=dpi)
