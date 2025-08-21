from settings import *

def plot_moments(setup, order, width_in_cm, height_in_cm):
    # Plot means
    moment_df = pd.read_csv(sim_results + f'moments/{setup}/sir_order{order}.csv')
    figsize = (width_in_cm/2.54, height_in_cm/2.54)
    # expected values
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(moment_df[(~moment_df.muS.isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df.muS.isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].muS, color=colors['S'], label='S')
    ax.plot(moment_df[~moment_df.muS.isna() & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df.muS.isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].muI, color=colors['I'], label='I')
    ax.plot(moment_df[~moment_df.muS.isna() & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df.muS.isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].muR, color=colors['R'], label='R')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Mean')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/{setup}/expected_value_sir{order}.png', dpi=dpi)

    if(order > 2):
        # Variances
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(moment_df[(~moment_df['M200'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M200'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)]['M200'], color=colors['S'], label='S')
        ax.plot(moment_df[(~moment_df['M020'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M020'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)]['M020'], color=colors['I'], label='I')
        ax.plot(moment_df[(~moment_df['M002'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M002'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)]['M002'], color=colors['R'], label='R')
        ax.legend()
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Variance')
        plt.tight_layout()
        fig.savefig(save_folder + f'moments/{setup}/variance_sir{order}.png', dpi=dpi)

        # Covariances
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(moment_df[(~moment_df['M110'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M110'].isna() & ((moment_df.muS <= moment_df.muS.iloc[0]))) & (moment_df.muS >= -5)]['M110'], color=colors['S'], label='S,I')
        ax.plot(moment_df[(~moment_df['M101'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M101'].isna() & ((moment_df.muS <= moment_df.muS.iloc[0]))) & (moment_df.muS >= -5)]['M101'], color=colors['I'], label='S,R')
        ax.plot(moment_df[(~moment_df['M011'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M011'].isna() & ((moment_df.muS <= moment_df.muS.iloc[0]))) & (moment_df.muS >= -5)]['M011'], color=colors['R'], label='I,R')
        ax.legend()
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Covariance')
        plt.tight_layout()
        fig.savefig(save_folder + f'moments/{setup}/covariance_sir{order}.png', dpi=dpi)

def plot_errors(setup, orders, width_in_cm, height_in_cm):
    moment_df_smm = pd.read_csv(sim_results + f'SMM/{setup}/moments.csv')
    figsize = (width_in_cm/2.54, height_in_cm/2.54)
    # expected value
    fig, ax = plt.subplots(figsize=figsize)
    colors = ['orange', 'blue', 'red', 'green', 'gray', 'black', 'teal', 'brown']
    for index, order in enumerate(orders):
        # Read moment df
        moment_df = pd.read_csv(sim_results + f'moments/{setup}/sir_order{order}.csv')
        # Merge smm and moment df
        df = pd.merge(moment_df_smm, moment_df, on="Time", how="outer", suffixes=("_df1", "_df2"))
        # Sort values by time
        df = df.sort_values("Time").reset_index(drop=True)
        # Interpolate missing values
        df = df.interpolate(method="linear")
        df["diff_muI"] = abs(df["muI_df1"] - df["muI_df2"])
        ax.plot(df.Time, df.diff_muI, color=colors[index], label=f'{order}-order')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Error')
    fig.suptitle('mu_I')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/{setup}/expected_value_error.png', dpi=dpi)
    
    # Variance
    fig, ax = plt.subplots(figsize=figsize)
    for index, order in enumerate(orders):
        if order == 2:
            continue
        # Read moment df
        moment_df = pd.read_csv(sim_results + f'moments/{setup}/sir_order{order}.csv')
        # Merge smm and moment df
        df = pd.merge(moment_df_smm, moment_df, on="Time", how="outer", suffixes=("_df1", "_df2"))
        # Sort values by time
        df = df.sort_values("Time").reset_index(drop=True)
        # Interpolate missing values
        df = df.interpolate(method="linear")
        df["diff_M020"] = abs(df["M020_df1"] - df["M020_df2"])
        ax.plot(df.Time, df["diff_M020"], color=colors[index], label=f'{order}-order')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Error')
    fig.suptitle('Var(I)')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/{setup}/var_error.png', dpi=dpi)
    
def plot_mean_var_infected(I0_list, width_in_cm, height_in_cm):
    figsize = (width_in_cm/2.54, height_in_cm/2.54)
    # expected values
    fig, ax = plt.subplots(figsize=figsize)
    for index, I0 in enumerate(I0_list):
        moment_df = pd.read_csv(sim_results + f'moments/I0={I0}/sir_order3.csv')
        ax.plot(moment_df[~moment_df.muS.isna() & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df.muS.isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].muI, color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('mu_I')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/muI_sir.png', dpi=dpi)
    
    # variance I
    fig, ax = plt.subplots(figsize=figsize)
    for index, I0 in enumerate(I0_list):
        moment_df = pd.read_csv(sim_results + f'moments/I0={I0}/sir_order3.csv')
        ax.plot(moment_df[(~moment_df['M020'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df['M020'].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)]['M020'], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Var(I)')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/VarI_sir.png', dpi=dpi)
    
    # variance S
    fig, ax = plt.subplots(figsize=figsize)
    col = 'M200'
    for index, I0 in enumerate(I0_list):
        moment_df = pd.read_csv(sim_results + f'moments/I0={I0}/sir_order3.csv')
        ax.plot(moment_df[(~moment_df[col].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df[col].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)][col], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Var(S)')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/VarS_sir.png', dpi=dpi)
    
    # Covariance S I
    fig, ax = plt.subplots(figsize=figsize)
    col = 'M110'
    for index, I0 in enumerate(I0_list):
        moment_df = pd.read_csv(sim_results + f'moments/I0={I0}/sir_order3.csv')
        ax.plot(moment_df[(~moment_df[col].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)].Time, moment_df[(~moment_df[col].isna()) & ((moment_df.muS <= moment_df.muS.iloc[0])) & (moment_df.muS >= -5)][col], color=other_colors[index], label=f'I0={I0}')
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Covar(S,I)')
    plt.tight_layout()
    fig.savefig(save_folder + f'moments/CovarSI_sir.png', dpi=dpi)
    
width_in_cm = 15
height_in_cm = 10
orders = [3]
I0_list = [1] #2, 3, 4, 5, 10, 50, 100

for I0 in I0_list:
    setup = f'I0={I0}RKStep'  # Initially infected individuals
    for order in orders:
        plot_moments(setup, order, width_in_cm, height_in_cm)

    plot_errors(setup, orders, width_in_cm, height_in_cm)
    
#plot_mean_var_infected(I0_list, width_in_cm, height_in_cm)
