from settings import *

num_sims = 500
dfs_smm = []
for i in range(num_sims):
    result_file = sim_results + f'SMM/{i}_comps.csv'
    df = pd.read_csv(result_file)
    dfs_smm.append(df)
    
# Assumes all have the same time points in the same order
time_smm = dfs_smm[0]['Time']
comps_smm = list(dfs_smm[0].columns)
comps_smm.remove('Time')

width_in_cm = 7
height_in_cm = 5.5

figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)
for c in comps_smm:
    matrix = np.array([df[c].values for df in dfs_smm])
    mean_val = np.mean(matrix, axis=0)
    p5 = np.percentile(matrix, 5, axis=0)
    p95 = np.percentile(matrix, 95, axis=0)
    ax.plot(time_smm, mean_val, color=colors[c], label=c)
    ax.fill_between(time_smm, p5, p95, color=colors[c], alpha=0.2, lw=0)
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + 'smm_sir.png', dpi=dpi)

dfs_abm = []
for i in range(num_sims):
    result_file = sim_results + f'ABM/{i}_comps.csv'
    df= pd.read_csv(result_file)
    df = df.drop(columns=['I_NS', 'I_Sy', 'I_Sev', 'I_Crit', 'D'])
    df = df.rename(columns={'E': 'I'})
    dfs_abm.append(df)
    
# Assumes all have the same time points in the same order
time_abm = dfs_abm[0]['Time']
comps = list(dfs_abm[0].columns)
comps.remove('Time')

figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)
for c in comps_smm:
    matrix_smm = np.array([df[c].values for df in dfs_smm])
    mean_val_smm = np.mean(matrix_smm, axis=0)
    p5_smm = np.percentile(matrix_smm, 5, axis=0)
    p95_smm = np.percentile(matrix_smm, 95, axis=0)
    matrix_abm = np.array([df[c].values for df in dfs_abm])
    mean_val_abm = np.mean(matrix_abm, axis=0)
    p5_abm = np.percentile(matrix_abm, 5, axis=0)
    p95_abm = np.percentile(matrix_abm, 95, axis=0)
    ax.plot(time_smm, mean_val_smm - mean_val_abm, color=colors[c], label=c)
    ax.fill_between(time_smm, p5_smm - p5_abm, p95_smm - p95_abm, color=colors[c], alpha=0.2, lw=0)
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + 'err_smm_abm_sir.png', dpi=dpi)

# Plot only I compartment for both models
figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)

matrix_smm = np.array([df['I'].values for df in dfs_smm])
mean_val_smm = np.mean(matrix_smm, axis=0)
p5_smm = np.percentile(matrix_smm, 5, axis=0)
p95_smm = np.percentile(matrix_smm, 95, axis=0)
matrix_abm = np.array([df['I'].values for df in dfs_abm])
mean_val_abm = np.mean(matrix_abm, axis=0)
p5_abm = np.percentile(matrix_abm, 5, axis=0)
p95_abm = np.percentile(matrix_abm, 95, axis=0)
ax.plot(time_smm, mean_val_smm, color=colors['S'], label='SMM')
ax.plot(time_abm, mean_val_abm, color=colors['I'], label='ABM')
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + 'smm_abm_I.png', dpi=dpi)
