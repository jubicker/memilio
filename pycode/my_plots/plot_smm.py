from settings import *

def plot_sim_outputs(setup, width_in_cm, height_in_cm, num_sims):
    # Plot simulation results (mean and 90 quantile)
    dfs_smm = []
    for i in range(num_sims):
        result_file = sim_results + f'SMM/{setup}/{i}_comps.csv'
        df = pd.read_csv(result_file)
        dfs_smm.append(df)
        
    # Assumes all have the same time points in the same order
    time_smm = dfs_smm[0]['Time']
    comps_smm = list(dfs_smm[0].columns)
    comps_smm.remove('Time')

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
    fig.savefig(save_folder + f'SMM/{setup}/smm_sir.png', dpi=dpi)

    # Plot means
    moment_df = pd.read_csv(sim_results + f'SMM/{setup}/moments.csv')
    figsize = (width_in_cm/2.54, height_in_cm/2.54)
    # expected values
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(moment_df.Time, moment_df.muS, color=colors['S'], label='S')
    ax.plot(moment_df.Time, moment_df.muI, color=colors['I'], label='I')
    ax.plot(moment_df.Time, moment_df.muR, color=colors['R'], label='R')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Mean')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/{setup}/mean_smm_sir.png', dpi=dpi)

    # Variances
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(moment_df.Time, moment_df['M200'], color=colors['S'], label='S')
    ax.plot(moment_df.Time, moment_df['M020'], color=colors['I'], label='I')
    ax.plot(moment_df.Time, moment_df['M002'], color=colors['R'], label='R')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Variance')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/{setup}/variance_smm_sir.png', dpi=dpi)

    # Covariances
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(moment_df.Time, moment_df['M110'], color=colors['S'], label='S,I')
    ax.plot(moment_df.Time, moment_df['M101'], color=colors['I'], label='S,R')
    ax.plot(moment_df.Time, moment_df['M011'], color=colors['R'], label='I,R')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Covariance')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/{setup}/covariance_smm_sir.png', dpi=dpi)

def plot_mean_var_infected(I0_list, width_in_cm, height_in_cm):
    figsize = (width_in_cm/2.54, height_in_cm/2.54)
    fig, ax = plt.subplots(figsize=figsize)
    # Means
    for index, I0 in enumerate(I0_list):
        moment_df = pd.read_csv(sim_results + f'SMM/I0={I0}/moments.csv')
        ax.plot(moment_df.Time, moment_df.muI, color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Mean Infected')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/mean_I.png', dpi=dpi)
    
    fig, ax = plt.subplots(figsize=figsize)
    # Variance I
    for index, I0 in enumerate(I0_list):
        # Plot means
        moment_df = pd.read_csv(sim_results + f'SMM/I0={I0}/moments.csv')
        # expected values
        ax.plot(moment_df.Time, moment_df['M020'], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Variance Infected')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/variance_I.png', dpi=dpi)
    
    fig, ax = plt.subplots(figsize=figsize)
    # Variance S
    for index, I0 in enumerate(I0_list):
        # Plot means
        moment_df = pd.read_csv(sim_results + f'SMM/I0={I0}/moments.csv')
        # expected values
        ax.plot(moment_df.Time, moment_df['M200'], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Variance Susceptible')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/variance_S.png', dpi=dpi)
    
    fig, ax = plt.subplots(figsize=figsize)
    # Covariance SI
    for index, I0 in enumerate(I0_list):
        # Plot means
        moment_df = pd.read_csv(sim_results + f'SMM/I0={I0}/moments.csv')
        # expected values
        ax.plot(moment_df.Time, moment_df['M110'], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Cov(S,I)')
    plt.tight_layout()
    fig.savefig(save_folder + f'SMM/Covariance_S_I.png', dpi=dpi)

width_in_cm = 15
height_in_cm = 10
num_sims = 1000

# for I0 in [100]:
#     setup = f'I0={I0}'  # Initially infected individuals
#     plot_sim_outputs(setup, width_in_cm, height_in_cm, num_sims)

plot_mean_var_infected([1,2,3,4,5,10,50,100], width_in_cm, height_in_cm)
