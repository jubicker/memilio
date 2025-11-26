from settings import *
import os

""" Plot all runs in one figure each for S, I, R """
def plot_all_runs(num_runs, result_dir, save_dir):
    figsize = (6,4)
    figS, axS = plt.subplots(figsize=figsize)
    figI, axI = plt.subplots(figsize=figsize)
    figR, axR = plt.subplots(figsize=figsize)
    for run in range(num_runs):
        data = pd.read_csv(result_dir + f"{run}_comps.csv")
        time = data["Time"]
        S = data.iloc[:,1]
        I = data.iloc[:,2]
        R = data.iloc[:,3]
        axS.plot(time, S)
        axI.plot(time, I)
        axR.plot(time, R)

    axS.set_xlabel("Time [days]")
    axI.set_xlabel("Time [days]")
    axR.set_xlabel("Time [days]")

    axS.set_ylabel("Susceptible [#]")
    axI.set_ylabel("Infected [#]")
    axR.set_ylabel("Recovered [#]")

    figS.savefig(save_dir + "all_runs_S.png", dpi=dpi)
    figI.savefig(save_dir + "all_runs_I.png", dpi=dpi)
    figR.savefig(save_dir + "all_runs_R.png", dpi=dpi)
    
    plt.tight_layout()

    plt.close(figS)
    plt.close(figI)
    plt.close(figR)

""" Plot percentiles (e.g. 5th and 95th) over all runs in one figure each for S, I, R, and all together """
def plot_percentiles(num_runs, result_dir, save_dir, percentiles):
    figsize = (6, 4)

    all_S, all_I, all_R = [], [], []
    time = []

    for run in range(num_runs):
        print(run)
        data = pd.read_csv(os.path.join(result_dir, f"{run}_comps.csv"))
        if len(time) == 0:
            time = data["Time"].values
        all_S.append(data.iloc[:, 1].values)
        all_I.append(data.iloc[:, 2].values)
        all_R.append(data.iloc[:, 3].values)

    S_arr = np.array(all_S)
    I_arr = np.array(all_I)
    R_arr = np.array(all_R)

    # Compute percentiles
    S_p1, S_p2 = np.percentile(S_arr, percentiles, axis=0)
    I_p1, I_p2 = np.percentile(I_arr, percentiles, axis=0)
    R_p1, R_p2 = np.percentile(R_arr, percentiles, axis=0)
    
    # Susceptible
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), S_p1, S_p2, color=compartment_colors["S"], alpha=0.3)
    ax.plot(np.array(time), np.mean(S_arr, axis=0), color=compartment_colors["S"])
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Susceptible [#]")
    fig.tight_layout()
    fig.savefig(save_dir + "percentiles_S.png", dpi=dpi)
    plt.close(fig)
    
    # Infected
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), I_p1, I_p2, color=compartment_colors["I"], alpha=0.3)
    ax.plot(np.array(time), np.mean(I_arr, axis=0), color=compartment_colors["I"])
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Infected [#]")
    fig.tight_layout()
    fig.savefig(save_dir + "percentiles_I.png", dpi=dpi)
    plt.close(fig)

    # Recovered
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), R_p1, R_p2, color=compartment_colors["R"], alpha=0.3)
    ax.plot(np.array(time), np.mean(R_arr, axis=0), color=compartment_colors["R"])
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Recovered [#]")
    fig.tight_layout()
    fig.savefig(save_dir + "percentiles_R.png", dpi=dpi)
    plt.close(fig)
    
    # All
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), S_p1, S_p2, color=compartment_colors["S"], alpha=0.3)
    ax.plot(np.array(time), np.mean(S_arr, axis=0), color=compartment_colors["S"], label="S")
    ax.fill_between(np.array(time), I_p1, I_p2, color=compartment_colors["I"], alpha=0.3)
    ax.plot(np.array(time), np.mean(I_arr, axis=0), color=compartment_colors["I"], label="I")
    ax.fill_between(np.array(time), R_p1, R_p2, color=compartment_colors["R"], alpha=0.3)
    ax.plot(np.array(time), np.mean(R_arr, axis=0), color=compartment_colors["R"], label="R")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Individuals [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "percentiles_all.png", dpi=dpi)
    plt.close(fig)

""" Plot mean and standard deviation for S, I, R """
def plot_mean_var(result_dir, save_dir):
    means = pd.read_csv(result_dir + f"means.csv")
    moments = pd.read_csv(result_dir + f"moments.csv")
    figsize = (6,4)
    time = means["Time"]
    total_pop = means['muS_r0'][0] + means['muI_r0'][0] + means['muR_r0'][0]
    
    # Susceptible
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muS_r0'], color=compartment_colors["S"], label = "Mean")
    ax.fill_between(time, means['muS_r0'] - np.sqrt(moments['M200']), means['muS_r0'] + np.sqrt(moments['M200']), color=compartment_colors["S"], alpha=0.5, label = "Stddev")
    ax.fill_between(time, means['muS_r0'] - 2*np.sqrt(moments['M200']), means['muS_r0'] + 2*np.sqrt(moments['M200']), color=compartment_colors["S"], alpha=0.3, label = "2Stddev")
    ax.fill_between(time, means['muS_r0'] - moments['M200'], means['muS_r0'] + moments['M200'], color=compartment_colors["S"], alpha=0.2, label = "Var")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Susceptible [#]")
    ax.set_ylim(-0.1*total_pop, total_pop)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "mean_var_S.png", dpi=dpi)
    plt.close(fig)
    #Infected
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muI_r0'], color=compartment_colors["I"], label = "Mean")
    ax.fill_between(time, means['muI_r0'] - np.sqrt(moments['M020']), means['muI_r0'] + np.sqrt(moments['M020']), color=compartment_colors["I"], alpha=0.5, label = "Stddev")
    ax.fill_between(time, means['muI_r0'] - 2*np.sqrt(moments['M020']), means['muI_r0'] + 2*np.sqrt(moments['M020']), color=compartment_colors["I"], alpha=0.3, label = "2Stddev")
    ax.fill_between(time, means['muI_r0'] - moments['M020'], means['muI_r0'] + moments['M020'], color=compartment_colors["I"], alpha=0.2, label = "Var")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Infected [#]")
    ax.set_ylim(-0.1*total_pop, total_pop)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "mean_var_I.png", dpi=dpi)
    plt.close(fig)
    #Recovered
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, means['muR_r0'], color=compartment_colors["R"], label = "Mean")
    ax.fill_between(time, means['muR_r0'] - np.sqrt(moments['M002']), means['muR_r0'] + np.sqrt(moments['M002']), color=compartment_colors["R"], alpha=0.5, label = "Stddev")
    ax.fill_between(time, means['muR_r0'] - 2*np.sqrt(moments['M002']), means['muR_r0'] + 2*np.sqrt(moments['M002']), color=compartment_colors["R"], alpha=0.3, label = "2Stddev")
    ax.fill_between(time, means['muR_r0'] - moments['M002'], means['muR_r0'] + moments['M002'], color=compartment_colors["R"], alpha=0.2, label = "Var")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Recovered [#]")
    ax.set_ylim(-0.1*total_pop, total_pop)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "mean_var_R.png", dpi=dpi)
    plt.close(fig)

""" Plot variances for S, I, R """
def plot_variances(result_dir, save_dir):
    moments = pd.read_csv(result_dir + f"moments.csv")
    figsize = (6,4)
    time = moments["Time"]
    
    values = {
        "Var(S)": moments['M200'],
        "Var(I)": moments['M020'],
        "Var(R)": moments['M002'],
    }
    
    fig, ax = plt.subplots(figsize=figsize)
    counter = 0
    for label, var in values.items():
        c = list(compartment_colors.keys())[counter]
        ax.plot(time, var, label=label, color=compartment_colors[c])
        fig1, ax1 = plt.subplots(figsize=figsize)
        ax1.plot(time, var, color=compartment_colors[c])
        ax1.set_xlabel("Time [days]")
        ax1.set_ylabel(label)
        fig1.tight_layout()
        fig1.savefig(save_dir + f"{label}.png", dpi=dpi)
        counter += 1
        plt.close(fig1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Variance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "all_variances.png", dpi=dpi)
    plt.close(fig)
    
""" Plot covariances for S, I, R """
def plot_covariances(result_dir, save_dir):
    moments = pd.read_csv(result_dir + f"moments.csv")
    figsize = (6,4)
    time = moments["Time"]
    
    values = {
        "Cov(S,I)": moments['M110'],
        "Cov(S,R)": moments['M101'],
        "Cov(I,R)": moments['M011'],
    }

    color_dict = {
        "Cov(S,I)": colors["dark blue"],
        "Cov(S,R)": colors["middle green"],
        "Cov(I,R)": colors["orange"],
    }

    fig, ax = plt.subplots(figsize=figsize)
    for label, var in values.items():
        ax.plot(time, var, label=label, color=color_dict[label])
        fig1, ax1 = plt.subplots(figsize=figsize)
        ax1.plot(time, var, color=color_dict[label])
        ax1.set_xlabel("Time [days]")
        ax1.set_ylabel(label)
        fig1.tight_layout()
        fig1.savefig(save_dir + f"{label}.png", dpi=dpi)
        plt.close(fig1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Covariance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "all_covariances.png", dpi=dpi)
    plt.close(fig)

""" Plot all moments of order three """
def plot_third_order_moments(result_dir, save_dir):
    moments = pd.read_csv(result_dir + f"moments.csv")
    figsize = (6,4)
    time = moments["Time"]
    
    values = {
        "M111": moments['M111'],
        "M210": moments['M210'],
        "M120": moments['M120'],
        "M201": moments['M201'],
        "M102": moments['M102'],
        "M021": moments['M021'],
        "M012": moments['M012'],
        "M300": moments['M300'],
        "M030": moments['M030'],
        "M003": moments['M003'],
    }

    color_dict = {
        "M111": colors["dark blue"],
        "M210": colors["middle green"],
        "M120": colors["orange"],
        "M201": colors["teal"],
        "M102": colors["purple"],
        "M021": colors["rose"],
        "M012": colors["dark red"],
        "M300": colors["brown"],
        "M030": colors["middle grey"],
        "M003": colors["red"],
    }

    fig, ax = plt.subplots(figsize=figsize)
    for label, var in values.items():
        ax.plot(time, var, label=label, color=color_dict[label])
        fig1, ax1 = plt.subplots(figsize=figsize)
        ax1.plot(time, var, color=color_dict[label])
        ax1.set_xlabel("Time [days]")
        ax1.set_ylabel(label)
        fig1.tight_layout()
        fig1.savefig(save_dir + f"{label}.png", dpi=dpi)
        plt.close(fig1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Moment value")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "all_third_order_moments.png", dpi=dpi)
    plt.close(fig)

dir = "V:/bick_ju/TemporalHybrid/Hybrid1/"
config = "config_1r"
switch_value = "switch_value_1.000000"
num_runs = 10000

result_dir = dir + config + "/" + switch_value + "/"
save_dir = "H:/Documents/TemporalHybridModel/Hybrid1/" + config + "/" + switch_value + "/"
os.makedirs(save_dir, exist_ok=True)

# plot_all_runs(num_runs, result_dir, save_dir)
# plot_percentiles(num_runs, result_dir, save_dir, [5, 95])
plot_mean_var(result_dir, save_dir)
plot_variances(result_dir, save_dir)
plot_covariances(result_dir, save_dir)
# plot_third_order_moments(result_dir, save_dir)
