from settings import *
import matplotlib.pyplot as plt

# dictionary has as first keys num agents and as second keys num runs
smm_times_1_core_SIR = {
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
smm_times_1_core_SIRS = {
    1000: {
        1: [0.00235284, 0.00224043, 0.00227182, 0.00221876, 0.0022186, 0.00224294],
        1000: [1.61145, 1.61674, 1.6168, 1.61356, 1.62268, 1.62141],
    },
    10000: {
        1: [0.0137948, 0.0135915, 0.0133272, 0.0133507, 0.0130476, 0.0135499],
        1000: [12.8582, 12.7138, 12.7483, 12.7284, 12.8161, 12.7937],
    },
    100000: {
        1: [0.118815, 0.120006, 0.120212, 0.136644, 0.117747, 0.129988],
        1000: [118.07, 124.992, 118.793, 118.186, 117.771, 118.046],
    },
    1000000: {
        1: [1.05472, 1.10437, 1.11181, 1.1205, 1.10172, 1.07843],
        1000: [1098.87, 1103.58, 1109.01, 1097.62, 1098.83, 1102.71],
    },
    10000000: {
        1: [10.3642, 10.314, 10.5115, 10.4145, 10.3868, 10.3418],
        1000: [10337.4],
    },
}
# dictionary has as first keys num regions and as second keys num runs
smm_regions_1_core_SIR = {
    1: {
        1: [0.0181881, 0.0179115, 0.0175461],
        1000: [16.8318, 17.0659, 16.8218],
    },
    2: {
        1: [0.0758597, 0.0771048, 0.0762401],
        1000: [74.8697, 74.4533, 74.2241],
    },
    4: {
        1: [0.389559, 0.383875, 0.379288],
        1000: [374.46, 374.454, 374.281],
    },
    8: {
        1: [1.93076, 1.9472, 1.93119],
        1000: [1787.84, 1787.89, 1791.8],
    },
    16: {
        1: [11.3264, 11.3237, 11.3099],
        1000: [10468.9, 10537.5, 10453.5],
    },
}
smm_regions_1_core_SIRS = {
    1: {
        1: [0.117373, 0.12181, 0.118487],
        1000: [119.134, 118.22, 118.145],
    },
    2: {
        1: [0.480712, 0.485542, 0.477801],
        1000: [479.407, 480.724, 483.322],
    },
    4: {
        1: [2.40026, 2.38976, 2.35469],
        1000: [2368.63, 2400.13, 2380.06],
    },
    8: {
        1: [11.6073, 11.5991, 11.6333],
        1000: [11533.2, 11483.6, 11481.4],
    },
    16: {
        1: [67.234, 67.3805, 67.3359],
        1000: [66385.1],
    },
}

smm_regions_fixed_pop_w_spatial_SIR = {
    1: {
        1: [],
        1000: [],
    },
    2: {
        1: [],
        1000: [],
    },
    4: {
        1: [0.958945, 0.970116, 0.952182],
        1000: [949.825, 958.031, 968.853],
    },
    8: {
        1: [2.5274, 2.51801, 2.54943],
        1000: [2434.29, 2418.46, 2426.69],
    },
    16: {
        1: [8.32002, 8.32631, 8.43601],
        1000: [7525.44, 7519.02, 7557.18],
    },
}

smm_regions_fixed_pop_w_spatial_SIRS = {
    1: {
        1: [],
        1000: [],
    },
    2: {
        1: [],
        1000: [],
    },
    4: {
        1: [],
        1000: [],
    },
    8: {
        1: [14.6015, 14.9211, 14.7448],
        1000: [14580.8, 14511.5, 14463.3],
    },
    16: {
        1: [43.4253, 43.7765, 43.9178],
        1000: [43367.6],
    },
}

smm_regions_fixed_pop_wo_spatial_SIR = {
    1: {
        1: [],
        1000: [],
    },
    2: {
        1: [],
        1000: [],
    },
    4: {
        1: [0.504712, 0.495513, 0.490131],
        1000: [484.961, 492.647, 494.545],
    },
    8: {
        1: [1.02852, 1.01227, 1.03463],
        1000: [883.416, 890.864, 871.154],
    },
    16: {
        1: [2.58079, 2.59119, 2.58511],
        1000: [1690.03, 1731.48, 1697.54],
    },
}

smm_regions_fixed_pop_wo_spatial_SIRS = {
    1: {
        1: [],
        1000: [],
    },
    2: {
        1: [],
        1000: [],
    },
    4: {
        1: [],
        1000: [],
    },
    8: {
        1: [5.31619, 5.27769, 5.33338],
        1000: [5124.53, 5180.15, 5150.35],
    },
    16: {
        1: [10.3339, 10.2974, 10.2825],
        1000: [9374.44, 9545.95, 9360.39],
    },
}


# dictionary has as first key num agnets and as second key integrator settings
moment_times_1_core_SIR = {
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

moment_times_1_core_SIRS = {
    1000: {
        r"adaptive $\Delta t$": [0.000965175, 0.000912255, 0.000898445, 0.000910435, 0.000906565, 0.000897986],
        r"adaptive $\Delta t_{max}=0.1$": [0.0215145, 0.0213532, 0.021674, 0.0214449, 0.0216336, 0.0215983],
    },
    10000: {
        r"adaptive $\Delta t$": [0.00158296, 0.00113684, 0.00108351, 0.00111839, 0.00108809, 0.0010838],
        r"adaptive $\Delta t_{max}=0.1$": [0.0210945, 0.0215784, 0.0215448, 0.0216074, 0.0213278, 0.0215885],
    },
    100000: {
        r"adaptive $\Delta t$": [0.0013486, 0.00130609, 0.00132178, 0.001292, 0.00130638, 0.00130805],
        r"adaptive $\Delta t_{max}=0.1$": [0.0226892, 0.0219435, 0.0216905, 0.0219941, 0.0221015, 0.022194],
    },
    1000000: {
        r"adaptive $\Delta t$": [0.00147051, 0.00143203, 0.00143698, 0.00142414, 0.00143347, 0.00143784],
        r"adaptive $\Delta t_{max}=0.1$": [0.0219897, 0.0212882, 0.0215868, 0.0215646, 0.0214043, 0.0213693],
    },
    10000000: {
        r"adaptive $\Delta t$": [0.00220207, 0.00168009, 0.0016345, 0.00164243, 0.00160726, 0.00168393],
        r"adaptive $\Delta t_{max}=0.1$": [0.0216463, 0.0214257, 0.0215064, 0.0214426, 0.0214875, 0.0216528],
    }
}

# dictionary has as first key num regions and as second key integrator settings
moment_regions_1_core_SIR = {
    1: {
        r"adaptive $\Delta t$": [0.00175172, 0.00120175, 0.00116055],
        r"adaptive $\Delta t_{max}=0.1$": [0.0224906, 0.0219943, 0.0222127],
    },
    2: {
        r"adaptive $\Delta t$": [0.0109056, 0.0107577, 0.0107707],
        r"adaptive $\Delta t_{max}=0.1$": [0.186444, 0.185511, 0.186899],
    },
    4: {
        r"adaptive $\Delta t$": [0.127962, 0.129134, 0.128576],
        r"adaptive $\Delta t_{max}=0.1$": [2.15922, 2.17243, 2.17413],
    },
    8: {
        r"adaptive $\Delta t$": [2.45793, 2.45989, 2.45731],
        r"adaptive $\Delta t_{max}=0.1$": [42.9197, 42.9597, 42.9653],
    },
    16: {
        r"adaptive $\Delta t$": [29.897, 29.9512, 29.9179],
        r"adaptive $\Delta t_{max}=0.1$": [521.651, 519.436, 520.208],
    }
}

moment_regions_1_core_SIRS = {
    1: {
        r"adaptive $\Delta t$": [0.00175003,  0.00127667, 0.00128534],
        r"adaptive $\Delta t_{max}=0.1$": [0.0226239, 0.0221367, 0.022219],
    },
    2: {
        r"adaptive $\Delta t$": [0.0132247, 0.012701, 0.0127133],
        r"adaptive $\Delta t_{max}=0.1$": [0.18866, 0.187676, 0.187194],
    },
    4: {
        r"adaptive $\Delta t$": [0.14912, 0.146246, 0.147061],
        r"adaptive $\Delta t_{max}=0.1$": [2.17958, 2.19909, 2.19778],
    },
    8: {
        r"adaptive $\Delta t$": [2.83246, 2.8339, 2.83812],
        r"adaptive $\Delta t_{max}=0.1$": [43.0892, 43.0761, 43.0599],
    },
    16: {
        r"adaptive $\Delta t$": [34.3173, 34.3332, 34.2922],
        r"adaptive $\Delta t_{max}=0.1$": [518.441, 519.019, 517.825],
    }
}

moment_regions_fixed_pop_w_spatial_SIR = {
    1: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    2: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    4: {
        r"adaptive $\Delta t$": [0.132485, 0.132729, 0.131971],
        r"adaptive $\Delta t_{max}=0.1$": [2.1816, 2.1835, 2.16238],
    },
    8: {
        r"adaptive $\Delta t$": [2.37744, 2.3834, 2.38373],
        r"adaptive $\Delta t_{max}=0.1$": [40.6301, 40.6358, 40.5548],
    },
    16: {
        r"adaptive $\Delta t$": [34.6261, 34.6163, 34.6531],
        r"adaptive $\Delta t_{max}=0.1$": [582.48, 582.921, 582.94],
    }
}

moment_regions_fixed_pop_w_spatial_SIRS = {
    1: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    2: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    4: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    8: {
        r"adaptive $\Delta t$": [2.71482, 2.70971, 2.7169],
        r"adaptive $\Delta t_{max}=0.1$": [40.6376, 40.5786, 40.5486],
    },
    16: {
        r"adaptive $\Delta t$": [35.7, 35.7233, 35.666],
        r"adaptive $\Delta t_{max}=0.1$": [586.613, 586.784, 586.053],
    }
}

moment_regions_fixed_pop_wo_spatial_SIR = {
    1: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    2: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    4: {
        r"adaptive $\Delta t$": [0.116952, 0.116822, 0.116385],
        r"adaptive $\Delta t_{max}=0.1$": [2.20128, 2.2023, 2.20751],
    },
    8: {
        r"adaptive $\Delta t$": [2.01798, 2.0218, 2.01475],
        r"adaptive $\Delta t_{max}=0.1$": [40.5487, 40.5139, 40.5747],
    },
    16: {
        r"adaptive $\Delta t$": [29.0344, 29.0406, 29.0516],
        r"adaptive $\Delta t_{max}=0.1$": [581.621, 581.879, 581.342],
    }
}

moment_regions_fixed_pop_wo_spatial_SIRS = {
    1: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    2: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    4: {
        r"adaptive $\Delta t$": [],
        r"adaptive $\Delta t_{max}=0.1$": [],
    },
    8: {
        r"adaptive $\Delta t$": [2.35522, 2.35766, 2.35288],
        r"adaptive $\Delta t_{max}=0.1$": [40.9162, 40.9673, 40.9691],
    },
    16: {
        r"adaptive $\Delta t$": [30.619, 30.6228, 30.619],
        r"adaptive $\Delta t_{max}=0.1$": [581.406, 581.625, 581.449],
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
            label=r"MoM full adaptive $\Delta t$", color=colors["purple"], marker="o")
    ax.plot(list(ode_dict.keys()), adaptive_restricted,
            label=r"MoM adaptive $\Delta t_{max}=0.1$", color=colors["rose"], marker="o")
    run1 = []
    runs1000 = []
    for x_value in stoch_dict.keys():
        run1.append(np.mean(stoch_dict[x_value][1]))
        runs1000.append(
            np.mean(stoch_dict[x_value][1000]))
    ax.plot(list(stoch_dict.keys()), run1,
            label=r"Stochastic $n_{sims}=1$", color=colors["dark teal"], marker="^")
    ax.plot(list(stoch_dict.keys()), runs1000,
            label=r"Stochastic $n_{sims}=1000$", color=colors["teal"], marker="^")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xticks(list(ode_dict.keys()))
    ax.set_xticklabels([rf"$10^{{{label}}}$" for label in x_labels])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Population size [#]")
    ax.set_ylabel("Runtime [s]")
    fig.subplots_adjust(left=0.16, bottom=0.2, top=0.98, right=0.98)
    fig.savefig(save_dir + "scaling.png", dpi=dpi)
    plt.close(fig)
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend.png", dpi=dpi)


def plot_region_scaling(ode_dict, stoch_dict, save_dir, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    adaptive_full = []
    adaptive_restricted = []
    x_labels = []
    for x_value in ode_dict.keys():
        x_labels.append(str(int(np.log2(x_value))))
        adaptive_full.append(
            np.mean(ode_dict[x_value][r"adaptive $\Delta t$"]))
        adaptive_restricted.append(
            np.mean(ode_dict[x_value][r"adaptive $\Delta t_{max}=0.1$"]))
    ax.plot(list(ode_dict.keys()), adaptive_full,
            label=r"MoM full adaptive $\Delta t$", color=colors["purple"], marker="o")
    ax.plot(list(ode_dict.keys()), adaptive_restricted,
            label=r"MoM adaptive $\Delta t_{max}=0.1$", color=colors["rose"], marker="o")
    run1 = []
    runs1000 = []
    for x_value in stoch_dict.keys():
        run1.append(np.mean(stoch_dict[x_value][1]))
        runs1000.append(
            np.mean(stoch_dict[x_value][1000]))
    ax.plot(list(stoch_dict.keys()), run1,
            label=r"Stochastic $n_{sims}=1$", color=colors["dark teal"], marker="^")
    ax.plot(list(stoch_dict.keys()), runs1000,
            label=r"Stochastic $n_{sims}=1000$", color=colors["teal"], marker="^")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xticks(list(ode_dict.keys()))
    ax.set_xticklabels([rf"$2^{{{label}}}$" for label in x_labels])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Regions [#]")
    ax.set_ylabel("Runtime [s]")
    fig.subplots_adjust(left=0.16, bottom=0.2, top=0.98, right=0.98)
    fig.savefig(save_dir + "scaling.png", dpi=dpi)
    plt.close(fig)
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=(1.5*figsize[0], figsize[1]))
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend.png", dpi=dpi)


save_dir_pop = "/Users/julia/sim_outputs/output/pop_scaling/OneCoreSIRS/"
save_dir_regions = "/Users/julia/sim_outputs/output/region_scaling/OneCoreSIRS/"
# fig_size = (5, 4)
plot_scaling(ode_dict=moment_times_1_core_SIRS,
             stoch_dict=smm_times_1_core_SIRS, save_dir=save_dir_pop, figsize=(5, 3.5))
plot_region_scaling(ode_dict=moment_regions_1_core_SIRS,
                    stoch_dict=smm_regions_1_core_SIRS, save_dir=save_dir_regions, figsize=(5, 3.5))
