from settings import *
import os

def plot_expected_values(smm_res, save_dir, moment_res, init_times):
    means = pd.read_csv(smm_res + f"means.csv")
    figsize = (6,4)
    time = means["Time"]
    df_moments = {}
    for t in init_times:
        ts = pd.read_csv(moment_res + f"{t}/expected_values.csv")
        df_moments[t] = ts
    
    # Susceptible
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muS_r0'], color='black', label = "Simulation mean", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.C1, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Susceptible [#]")
    ax.set_ylim(bottom=-10, top= max(means['muS_r0']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "S.png", dpi=dpi)
    plt.close(fig)
    
    # Infected
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muI_r0'], color='black', label = "Simulation mean", linestyle = "dashed")
    max_I = max(means['muI_r0'])
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        max_I = max(max_I, max(df_clean.C2))
        ax.plot(df_clean.Time, df_clean.C2, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Infected [#]")
    ax.set_ylim(bottom=-10, top= max_I * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "I.png", dpi=dpi)
    plt.close(fig)
    
    # Recovered
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muR_r0'], color='black', label = "Simulation mean", linestyle = "dashed")
    max_R = max(means['muR_r0'])
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        max_R = max(max_R, max(df_clean.C3))
        ax.plot(df_clean.Time, df_clean.C3, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Recovered [#]")
    ax.set_ylim(bottom=-10, top= max_R * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "R.png", dpi=dpi)
    plt.close(fig)
    
def plot_variances(smm_res, save_dir, moment_res, init_times):
    means = pd.read_csv(smm_res + f"moments.csv")
    figsize = (6,4)
    time = means["Time"]
    df_moments = {}
    for t in init_times:
        ts = pd.read_csv(moment_res + f"{t}/moments.csv")
        df_moments[t] = ts
    
    # Susceptible
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M200'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M200, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Var(S)")
    ax.set_ylim(bottom=-10, top= max(means['M200']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Var_S.png", dpi=dpi)
    plt.close(fig)
    
    # Infected
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M020'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M020, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Var(I)")
    ax.set_ylim(bottom=-10, top= max(means['M020']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Var_I.png", dpi=dpi)
    plt.close(fig)
    
    # Recovered
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M002'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M002, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Var(R)")
    ax.set_ylim(bottom=-10, top= max(means['M002']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Var_R.png", dpi=dpi)
    plt.close(fig)

def plot_covariances(smm_res, save_dir, moment_res, init_times):
    means = pd.read_csv(smm_res + f"moments.csv")
    figsize = (6,4)
    time = means["Time"]
    df_moments = {}
    for t in init_times:
        ts = pd.read_csv(moment_res + f"{t}/moments.csv")
        df_moments[t] = ts
    
    # S, I
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M110'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M110, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Cov(S,I)")
    ax.set_ylim(bottom=min(means['M110']) * 1.1, top= max(means['M110']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Cov_SI.png", dpi=dpi)
    plt.close(fig)
    
    # I, R
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M011'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M011, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Cov(I,R)")
    ax.set_ylim(bottom=min(means['M011']) * 1.1, top= max(means['M011']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Cov_IR.png", dpi=dpi)
    plt.close(fig)
    
    # S, R
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['M101'], color='black', label = "Simulation", linestyle = "dashed")
    for t in init_times:
        ts = df_moments[t]
        df_clean = ts.dropna()
        ax.plot(df_clean.Time, df_clean.M101, label=f"Moments t={t}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Cov(S,R)")
    ax.set_ylim(bottom=min(means['M101']) * 1.1, top= max(means['M101']) * 1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "Cov_SR.png", dpi=dpi)
    plt.close(fig)

dir_moments = "V:/bick_ju/TemporalHybrid/Moments/"
dir_smm = "V:/bick_ju/TemporalHybrid/SMM/"
config = "config_1r"
closure_order = "15"
init_time = ["0.000000", "10.000000", "15.000000", "20.000000", "25.000000"]

result_dir_moments = dir_moments + config + "/closure_order_" + closure_order + "/"
result_dir_smm = dir_smm + config + "/"
save_dir = "H:/Documents/TemporalHybridModel/Moments/" + config + "/closure_order_" + closure_order + "/"
os.makedirs(save_dir, exist_ok=True)

plot_expected_values(result_dir_smm, save_dir, result_dir_moments, init_time)
plot_variances(result_dir_smm, save_dir, result_dir_moments, init_time)
plot_covariances(result_dir_smm, save_dir, result_dir_moments, init_time)
