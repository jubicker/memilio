from settings import *

num_sims = 500
stoch_model = 'ABM'

dfs_ode = []
result_file = sim_results + f'ODE/sir.csv'
df = pd.read_csv(result_file)
dfs_ode.append(df)
    
# Assumes all have the same time points in the same order
time_ode = dfs_ode[0]['Time']
comps_ode = list(dfs_ode[0].columns)
comps_ode.remove('Time')

dfs_abm = []
for i in range(num_sims):
    result_file = sim_results + f'{stoch_model}/{i}_comps.csv'
    df= pd.read_csv(result_file)
    if(len(df.columns) > 4):
        df = df.drop(columns=['I_NS', 'I_Sy', 'I_Sev', 'I_Crit', 'D'])
        df = df.rename(columns={'E': 'I'})
    dfs_abm.append(df)
    
# Assumes all have the same time points in the same order
time_abm = dfs_abm[0]['Time']
comps = list(dfs_abm[0].columns)
comps.remove('Time')

width_in_cm = 7
height_in_cm = 5.5
figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)
for c in comps:
    matrix_ode = np.array([df[c].values for df in dfs_ode])
    mean_val_ode = np.mean(matrix_ode, axis=0)
    p5_ode = np.percentile(matrix_ode, 5, axis=0)
    p95_ode = np.percentile(matrix_ode, 95, axis=0)
    matrix_abm = np.array([df[c].values for df in dfs_abm])
    mean_val_abm = np.mean(matrix_abm, axis=0)
    p5_abm = np.percentile(matrix_abm, 5, axis=0)
    p95_abm = np.percentile(matrix_abm, 95, axis=0)
    ax.plot(time_ode, mean_val_ode - mean_val_abm, color=colors[c], label=c)
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + f'err_ODE_{stoch_model}_sir.png', dpi=dpi)

# Plot only I compartment for both models
figsize = (width_in_cm/2.54, height_in_cm/2.54)
fig, ax = plt.subplots(figsize=figsize)

matrix_smm = np.array([df['I'].values for df in dfs_ode])
mean_val_smm = np.mean(matrix_smm, axis=0)
p5_smm = np.percentile(matrix_smm, 5, axis=0)
p95_smm = np.percentile(matrix_smm, 95, axis=0)
matrix_abm = np.array([df['I'].values for df in dfs_abm])
mean_val_abm = np.mean(matrix_abm, axis=0)
p5_abm = np.percentile(matrix_abm, 5, axis=0)
p95_abm = np.percentile(matrix_abm, 95, axis=0)
ax.plot(time_ode, mean_val_smm, color=colors['S'], label='ODE')
ax.plot(time_abm, mean_val_abm, color=colors['I'], label=f'{stoch_model}')
ax.legend()
ax.set_xlabel('Time [days]')
ax.set_ylabel('Individuals [#]')
plt.tight_layout()
fig.savefig(save_folder + f'ODE_{stoch_model}_I.png', dpi=dpi)
