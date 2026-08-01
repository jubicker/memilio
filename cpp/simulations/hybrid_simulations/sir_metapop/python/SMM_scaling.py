from settings import *
import matplotlib.pyplot as plt

# dictionary has as first keys num agents and as second keys num runs
smm_times_1_core = {
    1000: {
        1: [0.00422288, 0.00415508, 0.00418733, 0.00436674, 0.00417868, 0.0042291],
        100: [0.0387601,  0.0387705, 0.0387321, 0.0385351, 0.038685, 0.0387536],
        1000: [0.367823, 0.370314, 0.368271, 0.367395, 0.367811, 0.372346],
    },
    10000: {
        1: [0.00605611, 0.00594334, 0.00598413, 0.00605261, 0.00598965, 0.00596847],
        100: [0.163953, 0.150785, 0.150253, 0.151814, 0.151254, 0.155072],
        1000: [1.54086, 1.54195, 1.53816, 1.52513, 1.5262, 1.50246],
    },
    100000: {
        1: [0.0214927, 0.0215178, 0.0212699, 0.0219186, 0.0213786, 0.0212221],
        100: [1.3095, 1.30536, 1.30127, 1.2876, 1.27572, 1.30072],
        1000: [12.7291, 13.0447, 13.0162, 12.9503, 13.061, 13.0828],
    },
    1000000: {
        1: [0.17165, 0.171462, 0.170684, 0.170641, 0.170532, 0.178345],
        100: [12.581, 12.5846, 12.8371, 12.833, 12.6985, 12.8543],
        1000: [126.864, 126.639, 126.423, 129.11, 125.907, 129.281],
    },
    10000000: {
        1: [1.65469, 1.65561, 1.65244, 1.65361, 1.64662, 1.65154],
        100: [126.617, 126.49, 126.221, 129.087, 126.445, 125.428],
        1000: [1268.32, 1262.8, 1284.84, 1269.51, 1251.19, 1258.19],
    },
}
# dictionary has as first key num agnets and as second key integrator settings
moment_times_1_core = {
    1000: {
        r"adaptive $\Delta t$": [0.00179063, 0.00124788, 0.00118606, 0.0011651, 0.00119818, 0.00116157],
        r"adaptive $\Delta t_{max}=0.1$": [0.0259234, 0.0253203, 0.0254147, 0.0253381, 0.0249482, 0.0253664],
        r"fixed $\Delta t=0.1$": [0.0256068, 0.0251347, 0.0253303, 0.0252013, 0.0250974, 0.0251725],
    },
    10000: {
        r"adaptive $\Delta t$": [0.00139621, 0.00133496, 0.00139941, 0.00136951, 0.00135934, 0.00138147],
        r"adaptive $\Delta t_{max}=0.1$": [0.0257225, 0.0251899, 0.025322, 0.0252465, 0.0251303, 0.0253503],
        r"fixed $\Delta t=0.1$": [0.0261432, 0.0256741, 0.025165, 0.0251148, 0.0256364, 0.0250841],
    },
    100000: {
        r"adaptive $\Delta t$": [0.00145307, 0.00144076, 0.00143607, 0.0014666, 0.00147857, 0.00144746],
        r"adaptive $\Delta t_{max}=0.1$": [0.0253991, 0.0256705, 0.0253495, 0.0255249, 0.0259488, 0.0255245],
        r"fixed $\Delta t=0.1$": [0.0255541, 0.0254212, 0.0252059, 0.0253196, 0.025387, 0.0253036],
    },
    1000000: {
        r"adaptive $\Delta t$": [0.00169478, 0.00162213, 0.00167822, 0.0016258, 0.00166331, 0.00159577],
        r"adaptive $\Delta t_{max}=0.1$": [0.0252094, 0.0253392, 0.0251693, 0.0253259, 0.0256441, 0.0257146],
        r"fixed $\Delta t=0.1$": [0.0253915, 0.0251108, 0.0252109, 0.0254144, 0.0251818, 0.0249921],
    },
    10000000: {
        r"adaptive $\Delta t$": [0.00183093, 0.00168246, 0.00176486, 0.00176509, 0.00164463, 0.00171075],
        r"adaptive $\Delta t_{max}=0.1$": [0.0254599, 0.0253249, 0.0252574, 0.0252438, 0.0254014, 0.0254253],
        r"fixed $\Delta t=0.1$": [0.0252072, 0.0252698, 0.0252897, 0.0251589, 0.0250423, 0.0253353],
    }
}


def plot_scaling(ode_dict, stoch_dict, save_dir, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    adaptive_full = []
    adaptive_restricted = []
    x_labels = []
    for x_value in ode_dict.keys():
        x_labels.append(str(int(np.log10(x_value))))
        adaptive_full.append(
            np.mean(ode_dict[x_value][r"adaptive $\Delta t$"]))
        adaptive_restricted.append(
            np.mean(ode_dict[x_value][r"adaptive $\Delta t_{max}=0.1$"]))
    ax.plot(list(ode_dict.keys()), adaptive_full,
            label=r"MoM full adaptive $\Delta t$", color=colors["dark blue"], marker="o")
    ax.plot(list(ode_dict.keys()), adaptive_restricted,
            label=r"MoM adaptive $\Delta t_{max}=0.1$", color=colors["middle blue"], marker="o")
    run1 = []
    runs1000 = []
    for x_value in stoch_dict.keys():
        run1.append(np.mean(stoch_dict[x_value][1]))
        runs1000.append(
            np.mean(stoch_dict[x_value][1000]))
    ax.plot(list(stoch_dict.keys()), run1,
            label=r"Stochastic $n_{sims}=1$", color=colors["dark green"], marker="^")
    ax.plot(list(stoch_dict.keys()), runs1000,
            label=r"Stochastic $n_{sims}=1000$", color=colors["middle green"], marker="^")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xticks(list(ode_dict.keys()))
    ax.set_xticklabels([rf"$10^{{{label}}}$" for label in x_labels])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Population size [#]")
    ax.set_ylabel("Runtime [s]")
    fig.subplots_adjust(left=0.15, bottom=0.18, top=0.98, right=0.98)
    fig.savefig(save_dir + "scaling.png", dpi=dpi)
    plt.close(fig)
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center')
    fig_leg.savefig(save_dir + "legend.png", dpi=dpi)


save_dir = "/Users/julia/sim_outputs/output/pop_scaling/OneCoreSIR/"
fig_size = (5, 3)
plot_scaling(ode_dict=moment_times_1_core,
             stoch_dict=smm_times_1_core, save_dir=save_dir, figsize=fig_size)
