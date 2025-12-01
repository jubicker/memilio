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
def plot_percentiles(result_dir, save_dir, percentiles, ode_dir = "", region = 0):
    figsize = (3.3, 2.5)

    all_S, all_I, all_R = [], [], []
    time = []
    
    # Read percentiles
    for p in percentiles:
        data = pd.read_csv(os.path.join(result_dir, f"{p}.csv"))
        all_S.append(data.iloc[:, region * 3 + 1].values)
        all_I.append(data.iloc[:, region * 3 + 2].values)
        all_R.append(data.iloc[:, region * 3 + 3].values)
        if len(time) == 0:
            time = data["Time"].values

    # Read mean
    mean_df = pd.read_csv(os.path.join(result_dir, f"means.csv"))
    mean_S = mean_df.iloc[:, region * 3 + 1].values
    mean_I = mean_df.iloc[:, region * 3 + 2].values
    mean_R = mean_df.iloc[:, region * 3 + 3].values
    
    # Susceptible
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), np.array(all_S[0]), np.array(all_S[1]), color=compartment_colors["S"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_S), color=compartment_colors["S"], label =  "Simulation")
    if ode_dir != "":
        # Read ode results
        ode_df = pd.read_csv(os.path.join(ode_dir, f"means.csv"))
        ode_S = ode_df.iloc[:, region * 3 + 1].values
        ax.plot(np.array(time), np.array(ode_S), color="black", linestyle="dashed", label = "ODE")      
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Susceptible [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"percentiles_S_{region}.png", dpi=dpi)
    plt.close(fig)
    
    # Infected
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), np.array(all_I[0]), np.array(all_I[1]), color=compartment_colors["I"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_I), color=compartment_colors["I"], label =  "Simulation")
    if ode_dir != "":
        # Read ode results
        ode_df = pd.read_csv(os.path.join(ode_dir, f"means.csv"))
        ode_I = ode_df.iloc[:, region * 3 + 2].values
        ax.plot(np.array(time), np.array(ode_I), color="black", linestyle="dashed", label = "ODE")      
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Infected [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"percentiles_I_{region}.png", dpi=dpi)
    plt.close(fig)

    # Recovered
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), np.array(all_R[0]), np.array(all_R[1]), color=compartment_colors["R"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_R), color=compartment_colors["R"], label =  "Simulation")
    if ode_dir != "":
        # Read ode results
        ode_df = pd.read_csv(os.path.join(ode_dir, f"means.csv"))
        ode_R = ode_df.iloc[:, region * 3 + 3].values
        ax.plot(np.array(time), np.array(ode_R), color="black", linestyle="dashed", label = "ODE")      
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Recovered [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"percentiles_R_{region}.png", dpi=dpi)
    plt.close(fig)
    
    # All
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(np.array(time), all_S[0], all_S[1], color=compartment_colors["S"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_S), color=compartment_colors["S"], label="S")
    ax.fill_between(np.array(time), all_I[0], all_I[1], color=compartment_colors["I"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_I), color=compartment_colors["I"], label="I")
    ax.fill_between(np.array(time), all_R[0], all_R[1], color=compartment_colors["R"], alpha=0.3)
    ax.plot(np.array(time), np.array(mean_R), color=compartment_colors["R"], label="R")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Individuals [#]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"percentiles_all_{region}.png", dpi=dpi)
    plt.close(fig)
    
def count_extinctions(num_runs, result_dir):
    num_extinctions = 0.
    for run in range(num_runs):
        print(run)
        data = pd.read_csv(os.path.join(result_dir, f"{run}_comps.csv"))
        if float(data.iloc[-1, 3]) < 200.:
            num_extinctions += 1
    print(f"Num extinctions: {num_extinctions}, Percentage: {num_extinctions/num_runs}")           

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
    ax.set_ylim(-0.1*total_pop, total_pop*1.1)
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
    ax.set_ylim(-0.1*total_pop, total_pop*1.1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "mean_var_R.png", dpi=dpi)
    plt.close(fig)

""" Plot variances for S, I, R """
def plot_variances(result_dir, save_dir, labels, colors):
    moment_list = []
    for res in result_dir:        
        moment_list.append(pd.read_csv(res + f"moments.csv"))
    figsize = (4, 3)
    time = moment_list[0]["Time"]
    
    values = {
        "Var(S)": [moments['M200'] for moments in moment_list],
        "Var(I)": [moments['M020'] for moments in moment_list],
        "Var(R)": [moments['M002'] for moments in moment_list],
    }
    
    # Plot variances
    fig, ax = plt.subplots(figsize=figsize)
    counter = 0
    for label, var in values.items():
        fig1, ax1 = plt.subplots(figsize=figsize)
        for config in range(len(moment_list)):
            ax.plot(time, var[config], label=label)
            ax1.plot(time, var[config], color=colors[config], label = labels[config])
            counter += 1
        ax1.set_xlabel("Time [days]")
        ax1.set_ylabel(label)
        ax1.legend()
        fig1.tight_layout()
        fig1.savefig(save_dir + f"{label}.png", dpi=dpi)
        plt.close(fig1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Variance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "all_variances.png", dpi=dpi)
    plt.close(fig)
    
    # Plot standard deviation
    values = {
        "Stddev(S)": [moments['M200'] for moments in moment_list],
        "Stddev(I)": [moments['M020'] for moments in moment_list],
        "Stddev(R)": [moments['M002'] for moments in moment_list],
    }
    
    fig, ax = plt.subplots(figsize=figsize)
    counter = 0
    for label, var in values.items():
        fig1, ax1 = plt.subplots(figsize=figsize)
        for config in range(len(moment_list)):
            ax.plot(time, np.sqrt(var[config]), label=label)
            ax1.plot(time, np.sqrt(var[config]), color=colors[config], label = labels[config])
            counter += 1
        ax1.set_xlabel("Time [days]")
        ax1.set_ylabel(label)
        ax1.legend()
        fig1.tight_layout()
        fig1.savefig(save_dir + f"{label}.png", dpi=dpi)
        plt.close(fig1)
        
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Stddev")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + "all_stddev.png", dpi=dpi)
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
config = ["config_2r1"]#["config1_1r_I0_1", "config1_1r_I0_2", "config1_1r_I0_10", "config1_1r_I0_100"]#
labels = ["I0=1", "I0=2", "I0=10", "I0=100"]
switch_value = "switch_value_1.000000"
num_runs = 10000

result_dir = []
for c in config:
    result_dir.append(dir + c + "/" + switch_value + "/")
save_dir = "H:/Documents/TemporalHybridModel/Hybrid1/" + config[0] + "/" + switch_value + "/"
#save_dir = "H:/Documents/TemporalHybridModel/Hybrid1/"
ode_dir = dir + config[0] + "/switch_value_0.000000/"
os.makedirs(save_dir, exist_ok=True)

# plot_all_runs(num_runs, result_dir, save_dir)
# plot_percentiles(num_runs, result_dir, save_dir, ["p05", "p95"], ode_dir)
# plot_mean_var(result_dir, save_dir)
# plot_variances(result_dir, save_dir, labels, [colors["rose"], colors["brown"], colors["dark blue"], colors["dark green"]])
# plot_covariances(result_dir, save_dir)
# plot_third_order_moments(result_dir, save_dir)

plot_percentiles(result_dir[0], save_dir, ["p05", "p95"], ode_dir, 0)
plot_percentiles(result_dir[0], save_dir, ["p05", "p95"], ode_dir, 1)

# count_extinctions(num_runs, result_dir[2])
