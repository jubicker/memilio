from settings import *
import matplotlib.pyplot as plt

pop_sizes = [1000, 10000, 100000, 1000000, 10000000]
region_size = [1, 2, 4, 8, 16]
linear_pop = [10**-3, 10**-2, 10**-1, 10**0, 10**1]
quadratic_pop = [10**-3, 10**-1, 10**1, 10**3, 10**5]
linear_regions = [10**-3, (2)*10**-3, (4)*10**-3, (8)*10**-3, (16)*10**-3]
quadratic_regions = [10**-3, (2**2)*10**-3, (4**2)
                     * 10**-3, (8**2)*10**-3, (16**2)*10**-3]
moment_regions = [10**-3, 2.62*10**-3, 3.03*2.62*10**-
                  3, 3.39*3.03*2.62*10**-3, 3.65*3.39*3.03*2.62*10**-3]
moment_quadratic_regions = [10**-3, (2**2*2.62)*10**-3, (4**2*3.03*2.62)
                            * 10**-3, (8**2*3.39*3.03*2.62)*10**-3, (16**2*3.65*3.39*3.03*2.62)*10**-3]

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
        1: [0.164895, 0.164689, 0.164669],
        1000: [162.921, 163.625, 163.891],
    },
    2: {
        1: [0.380913, 0.375082, 0.373749],
        1000: [379.39, 373.943, 374.27],
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
        1: [1.07152, 1.11267, 1.0895],
        1000: [1095.62, 1114.32, 1097.04],
    },
    2: {
        1: [2.36861, 2.35022, 2.31826],
        1000: [2469.5, 2376.01, 2359.96],
    },
    4: {
        1: [5.78239, 5.79902, 6.03084],
        1000: [5779.27, 5783.97, 5856.48],
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
        1: [0.164895, 0.164689, 0.164669],
        1000: [162.921, 163.625, 163.891],
    },
    2: {
        1: [0.289383, 0.284427, 0.288879],
        1000: [284.06, 283.448, 283.722],
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
        1: [1.07152, 1.11267, 1.0895],
        1000: [1095.62, 1114.32, 1097.04],
    },
    2: {
        1: [1.79085, 1.78082, 1.80709],
        1000: [1776.74, 1788.44, 1779.26],
    },
    4: {
        1: [3.12027, 3.06646, 3.04975],
        1000: [2965.94, 2956.7, 3003.17],
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
        r"adaptive $\Delta t$": [0.000912166, 0.000847095, 0.000849015],
        r"adaptive $\Delta t_{max}=0.1$": [0.0205028, 0.0197733, 0.0198291],
    },
    10000: {
        r"adaptive $\Delta t$": [0.00100269, 0.00098405, 0.000967711],
        r"adaptive $\Delta t_{max}=0.1$": [0.0197418, 0.0196858, 0.0196099],
    },
    100000: {
        r"adaptive $\Delta t$": [0.00109466, .00107103, 0.00106616],
        r"adaptive $\Delta t_{max}=0.1$": [0.0196336, 0.0196551, 0.0195718],
    },
    1000000: {
        r"adaptive $\Delta t$": [0.00118854, 0.00114375, 0.00114727],
        r"adaptive $\Delta t_{max}=0.1$": [0.0199979, 0.0196208, 0.0198075],
    },
    10000000: {
        r"adaptive $\Delta t$": [0.00122279, 0.00119657, 0.00120483],
        r"adaptive $\Delta t_{max}=0.1$": [0.0195449, 0.0197373, 0.0194327],
    }
}

moment_times_1_core_SIRS = {
    1000: {
        r"adaptive $\Delta t$": [0.000860244, 0.000818103, 0.000827784],
        r"adaptive $\Delta t_{max}=0.1$": [0.0197402, 0.0198888, 0.0197024],
    },
    10000: {
        r"adaptive $\Delta t$": [0.0010365, 0.000970914, 0.000972674],
        r"adaptive $\Delta t_{max}=0.1$": [0.0195802, 0.0198618, 0.0194906],
    },
    100000: {
        r"adaptive $\Delta t$": [0.00120201, 0.00117634, 0.0011618],
        r"adaptive $\Delta t_{max}=0.1$": [0.0197018, 0.0198889, 0.0197945],
    },
    1000000: {
        r"adaptive $\Delta t$": [0.0013366, 0.00133855, 0.00134141],
        r"adaptive $\Delta t_{max}=0.1$": [0.0195958, 0.0196706, 0.0197181],
    },
    10000000: {
        r"adaptive $\Delta t$": [0.00160484, 0.00154687, 0.00156641],
        r"adaptive $\Delta t_{max}=0.1$": [0.0197281, 0.0196126, 0.0198854],
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
        r"adaptive $\Delta t$": [0.00114108, 0.00111319, 0.00112762],
        r"adaptive $\Delta t_{max}=0.1$": [0.0193596, 0.0194053, 0.0193037],
    },
    2: {
        r"adaptive $\Delta t$": [0.0102234, 0.010259, 0.010418],
        r"adaptive $\Delta t_{max}=0.1$": [0.16247, 0.162177, 0.161401],
    },
    4: {
        r"adaptive $\Delta t$": [0.113578, 0.114378, 0.113785],
        r"adaptive $\Delta t_{max}=0.1$": [1.88329, 1.8878, 1.89014],
    },
    8: {
        r"adaptive $\Delta t$": [2.00702, 2.0075, 2.00724],
        r"adaptive $\Delta t_{max}=0.1$": [34.3889, 34.3955, 34.405],
    },
    16: {
        r"adaptive $\Delta t$": [23.3569, 23.3567, 23.3439],
        r"adaptive $\Delta t_{max}=0.1$": [392.654, 392.716, 392.66],
    }
}

moment_regions_fixed_pop_w_spatial_SIRS = {
    1: {
        r"adaptive $\Delta t$": [0.00182819, 0.00128907, 0.00127542],
        r"adaptive $\Delta t_{max}=0.1$": [0.0198133, 0.0192949, 0.0196479],
    },
    2: {
        r"adaptive $\Delta t$": [0.0115611, 0.011624, 0.0116741],
        r"adaptive $\Delta t_{max}=0.1$": [0.161299, 0.161371, 0.162185],
    },
    4: {
        r"adaptive $\Delta t$": [0.133818, 0.134561, 0.13405],
        r"adaptive $\Delta t_{max}=0.1$": [1.87947, 1.89004, 1.88599],
    },
    8: {
        r"adaptive $\Delta t$": [2.28469, 2.28418, 2.28524],
        r"adaptive $\Delta t_{max}=0.1$": [34.409, 34.407, 34.389],
    },
    16: {
        r"adaptive $\Delta t$": [24.0277, 23.9882, 24.0323],
        r"adaptive $\Delta t_{max}=0.1$": [392.247, 392.335, 392.901],
    }
}

moment_regions_fixed_pop_wo_spatial_SIR = {
    1: {
        r"adaptive $\Delta t$": [0.00172833, 0.00114751, 0.00113551],
        r"adaptive $\Delta t_{max}=0.1$": [0.0197243, 0.0195539, 0.0196301],
    },
    2: {
        r"adaptive $\Delta t$": [0.00642788, 0.00580122, 0.00574156],
        r"adaptive $\Delta t_{max}=0.1$": [0.108306, 0.107967, 0.108339],
    },
    4: {
        r"adaptive $\Delta t$": [0.0334687, 0.0336104, 0.0335257],
        r"adaptive $\Delta t_{max}=0.1$": [0.633772, 0.633355, 0.637646],
    },
    8: {
        r"adaptive $\Delta t$": [0.352347, 0.351245, 0.352093],
        r"adaptive $\Delta t_{max}=0.1$": [7.11487, 7.14482, 7.14497],
    },
    16: {
        r"adaptive $\Delta t$": [2.57063, 2.57961, 2.58441],
        r"adaptive $\Delta t_{max}=0.1$": [51.8698, 51.8843, 51.8754],
    }
}

moment_regions_fixed_pop_wo_spatial_SIRS = {
    1: {
        r"adaptive $\Delta t$": [0.00183568, 0.00132248, 0.00134133],
        r"adaptive $\Delta t_{max}=0.1$": [0.0196285, 0.0197345, 0.019624],
    },
    2: {
        r"adaptive $\Delta t$": [0.00705451, 0.00702802, 0.0067509],
        r"adaptive $\Delta t_{max}=0.1$": [0.109635, 0.109349, 0.109183],
    },
    4: {
        r"adaptive $\Delta t$": [0.0368772, 0.0369805, 0.0369522],
        r"adaptive $\Delta t_{max}=0.1$": [0.627141, 0.63235, 0.631196],
    },
    8: {
        r"adaptive $\Delta t$": [0.409035, 0.409538, 0.408569],
        r"adaptive $\Delta t_{max}=0.1$": [7.08212, 7.10519, 7.10707],
    },
    16: {
        r"adaptive $\Delta t$": [2.71083, 2.70973, 2.70624],
        r"adaptive $\Delta t_{max}=0.1$": [51.8427, 51.8316, 51.8781],
    }
}


def save_legends(ax, save_dir, figsize, complexity_file):
    # the complexity reference lines (O(...)) get their own legend png
    handles, labels = ax.get_legend_handles_labels()
    is_complexity = [r"\mathcal{O}" in label for label in labels]
    for file_name, select in [("legend.png", False), (complexity_file, True)]:
        sel_handles = [h for h, c in zip(
            handles, is_complexity) if c == select]
        sel_labels = [l for l, c in zip(labels, is_complexity) if c == select]
        fig_leg = plt.figure(figsize=(figsize[0]*3, figsize[1]))
        fig_leg.legend(sel_handles, sel_labels, loc='center',
                       ncol=len(sel_labels))
        fig_leg.savefig(save_dir + file_name, dpi=dpi)
        plt.close(fig_leg)


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

    run1 = []
    runs1000 = []
    for x_value in stoch_dict.keys():
        run1.append(np.mean(stoch_dict[x_value][1]))
        runs1000.append(
            np.mean(stoch_dict[x_value][1000]))

    # set y-ticks at every power of 10
    all_y = np.concatenate([np.array(adaptive_full),
                            np.array(adaptive_restricted),
                            np.array(run1),
                            np.array(runs1000)])
    all_y = all_y[all_y > 0]
    ax.set_yscale("log")
    ax.set_xscale("log")
    if all_y.size > 0:
        ymin_pow = int(np.floor(np.log10(all_y.min())))
        ymax_pow = int(np.ceil(np.log10(all_y.max())))
        yticks = [10.0 ** i for i in range(ymin_pow, ymax_pow + 1)]
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            [rf"$10^{{{i}}}$" for i in range(ymin_pow, ymax_pow + 1)])

    ax.plot(list(ode_dict.keys()), adaptive_full,
            label=r"MoM full adaptive $\Delta t$", color=colors["purple"], marker="o")
    ax.plot(list(ode_dict.keys()), adaptive_restricted,
            label=r"MoM adaptive $\Delta t_{max}=0.1$", color=colors["rose"], marker="o")
    ax.plot(list(stoch_dict.keys()), run1,
            label=r"Stochastic $n_{sims}=1$", color=colors["dark teal"], marker="^")
    ax.plot(list(stoch_dict.keys()), runs1000,
            label=r"Stochastic $n_{sims}=1000$", color=colors["teal"], marker="^")

    ax.plot(pop_sizes, linear_pop,
            color=colors["dark grey"], linestyle="dashed", alpha=0.8, label=r"$\mathcal{O}(N)$")
    ax.plot(pop_sizes, quadratic_pop,
            color=colors["dark grey"], linestyle="dotted", alpha=0.8, label=r"$\mathcal{O}(N^2)$")

    ax.set_ylim(10**(-3)-10, 10**(5.2))
    ax.set_xticks(list(ode_dict.keys()))
    ax.set_xticklabels([rf"$10^{{{label}}}$" for label in x_labels])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Population size [#]")
    ax.set_ylabel("Runtime [s]")
    fig.subplots_adjust(left=0.17, bottom=0.2, top=0.97, right=0.98)
    fig.savefig(save_dir + "scaling_p.png", dpi=dpi)
    plt.close(fig)
    save_legends(ax, save_dir, figsize, "legend_complexity_p.png")


def plot_region_scaling(ode_dict, stoch_dict, save_dir, figsize, x_postfix):
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
    run1 = []
    runs1000 = []
    for x_value in stoch_dict.keys():
        run1.append(np.mean(stoch_dict[x_value][1]))
        runs1000.append(
            np.mean(stoch_dict[x_value][1000]))
    ax.set_yscale("log")
    ax.set_xscale("log")

    # set y-ticks at every power of 10
    all_y = np.concatenate([np.array(adaptive_full),
                            np.array(adaptive_restricted),
                            np.array(run1),
                            np.array(runs1000)])
    all_y = all_y[all_y > 0]
    if all_y.size > 0:
        ymin_pow = int(np.floor(np.log10(all_y.min())))
        ymax_pow = int(np.ceil(np.log10(all_y.max())))
        yticks = [10.0 ** i for i in range(ymin_pow, ymax_pow + 3)]
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            [rf"$10^{{{i}}}$" for i in range(ymin_pow, ymax_pow + 3)])

    ax.plot(list(ode_dict.keys()), adaptive_full,
            label=r"MoM full adaptive $\Delta t$", color=colors["purple"], marker="o")
    ax.plot(list(ode_dict.keys()), adaptive_restricted,
            label=r"MoM adaptive $\Delta t_{max}=0.1$", color=colors["rose"], marker="o")
    ax.plot(list(stoch_dict.keys()), run1,
            label=r"Stochastic $n_{sims}=1$", color=colors["dark teal"], marker="^")
    ax.plot(list(stoch_dict.keys()), runs1000,
            label=r"Stochastic $n_{sims}=1000$", color=colors["teal"], marker="^")

    ax.plot(region_size, linear_regions,
            color=colors["dark grey"], linestyle=(0, (5, 10)), alpha=0.8, label=r"$\mathcal{O}(3n_R)$")
    ax.plot(region_size, quadratic_regions,
            color=colors["dark grey"], linestyle="dashed", alpha=0.8, label=r"$\mathcal{O}(3n_R^2)$")
    ax.plot(region_size, moment_regions,
            color=colors["dark grey"], linestyle="dotted", alpha=0.8, label=r"$\mathcal{O}(n_m+3n_R)$")
    ax.plot(region_size, moment_quadratic_regions,
            color=colors["dark grey"], linestyle="dashdot", alpha=0.8, label=r"$\mathcal{O}(3n_R^2(n_m+3n_R))$")

    ax.set_xticks(list(ode_dict.keys()))
    ax.set_xticklabels([rf"$2^{{{label}}}$" for label in x_labels])
    ax.grid(visible=True, color=colors["middle grey"],
            linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_xlabel(f"Regions ({x_postfix}) [#]")
    ax.set_ylabel("Runtime [s]")
    ax.set_ylim(10**(-3)-10, 10**(5.2))
    fig.subplots_adjust(left=0.17, bottom=0.2, top=0.97, right=0.98)
    fig.savefig(save_dir + "scaling_r.png", dpi=dpi)
    plt.close(fig)
    save_legends(ax, save_dir, figsize, "legend_complexity_r.png")


save_dir_pop = ""
save_dir_regions = ""
# fig_size = (5, 4)
plot_scaling(ode_dict=moment_times_1_core_SIR,
             stoch_dict=smm_times_1_core_SIR, save_dir=save_dir_pop, figsize=(5, 3.5))
plot_region_scaling(ode_dict=moment_regions_fixed_pop_w_spatial_SIR,
                    stoch_dict=smm_regions_fixed_pop_w_spatial_SIR, save_dir=save_dir_regions + "w_spatial", figsize=(5, 3.5), x_postfix=r"$\kappa_z^{(k,l)}>0$")
plot_region_scaling(ode_dict=moment_regions_fixed_pop_wo_spatial_SIR,
                    stoch_dict=smm_regions_fixed_pop_wo_spatial_SIR, save_dir=save_dir_regions + "wo_spatial", figsize=(5, 3.5), x_postfix=r"$\kappa_z^{(k,l)}=0$")
