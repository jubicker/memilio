import pandas as pd
import matplotlib.pyplot as plt
import os
from settings import *

def plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index, figsize, tmax, num_regions):
    fig_var, ax_var = plt.subplots(figsize=figsize)
    fig_std, ax_std = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        moments = pd.read_csv(dir + f"/{config}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        time = moments['Time']
        ax_var.plot(time, moments[col_name], label=config_labels[c], color = compartment_colors[comp][c])
        ax_std.plot(time, np.sqrt(moments[col_name]), label=config_labels[c], color = compartment_colors[comp][c])

    ax_var.set_xlabel("Time [days]")
    ax_var.set_ylabel(f"Var {compartment_names[comp]}")
    ax_var.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax_var.legend()
    fig_var.tight_layout()
    fig_var.savefig(save_dir + f"/all_variances_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_var)
    
    ax_std.set_xlabel("Time [days]")
    ax_std.set_ylabel(f"Stddev {compartment_names[comp]}")
    ax_std.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax_std.legend()
    fig_std.tight_layout()
    fig_std.savefig(save_dir + f"/all_standard_deviations_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_std)
    
def plot_variance_and_std_two_models(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, model1_name, model2_name, save_dir, comp_index, figsize, tmax, num_regions, colors_model1, colors_model2):
    fig_var, ax_var = plt.subplots(figsize=figsize)
    fig_std, ax_std = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        moments = pd.read_csv(dir_model1 + f"/{config}/{condition}/{closure_method}/{closure_order}/{init_time}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        time = moments['Time']
        ax_var.plot(time, moments[col_name], label=model1_name + " " + config_labels[c], color = colors_model1[c], linewidth=2)
        ax_std.plot(time, np.sqrt(moments[col_name]), label=model1_name + " " + config_labels[c], color = colors_model1[c])
        
        moments = pd.read_csv(dir_model2 + f"/{config}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        time = moments['Time']
        ax_var.plot(time, moments[col_name], label=model2_name + " " + config_labels[c], color = colors_model2[c], linestyle="dotted", linewidth=2)
        ax_std.plot(time, np.sqrt(moments[col_name]), label=model2_name + " " + config_labels[c], color = colors_model2[c], linestyle="dotted")

    ax_var.set_xlabel("Time [days]")
    ax_var.set_ylabel(r"$\sigma^2_I$")
    #ax_var.set_yscale("log")
    ax_var.ticklabel_format(axis='y', style='sci', scilimits=(11,11))
    # if(config_labels[0] > ""):
    #     ax_var.legend()
    fig_var.tight_layout()
    fig_var.savefig(save_dir + f"/all_variances_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_var)
    
    ax_std.set_xlabel("Time [days]")
    ax_std.set_ylabel(r"$\sigma_I$")
    ax_std.ticklabel_format(axis='y', style='sci', scilimits=(11,11))
    # if(config_labels[0] > ""):
    #     ax_std.legend()
    fig_std.tight_layout()
    fig_std.savefig(save_dir + f"/all_standard_deviations_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig_std)

def plot_means_two_models(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, model1_name, model2_name, save_dir, comp_index, figsize, tmax, num_regions, colors_model1, colors_model2):
    fig, ax = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    for c, config in enumerate(configs):        
        means = pd.read_csv(dir_model1 + f"/{config}/{condition}/{closure_method}/{closure_order}/{init_time}/means.csv")
        means = means[means['Time'] <= tmax]
        time = means['Time']
        ax.plot(time, means.iloc[:, 1 + comp_index], label=model1_name + " " + config_labels[c], color = colors_model1[c], linewidth=2)
        
        means = pd.read_csv(dir_model2 + f"/{config}/means.csv")
        print(config, np.max(means.iloc[:, 1 + comp_index])/10000000)
        means = means[means['Time'] <= tmax]
        time = means['Time']
        ax.plot(time, means.iloc[:, 1 + comp_index], label=model2_name + " " + config_labels[c], color = colors_model2[c], linestyle="dotted", linewidth=2)

    ax.set_xlabel("Time [days]")
    ax.set_ylabel(r"$\mu_I$")
    #ax.set_yscale('log')
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    # if(config_labels[0] > ""):
    #     ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"/all_means_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/all_meanslegend.png", dpi=dpi, bbox_inches='tight')
    plt.close(fig_leg)
    plt.close(fig)
    
def plot_means(dir, configs, config_labels, save_dir, comp_index, figsize, tmax, num_regions):
    fig, ax = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    for c, config in enumerate(configs):        
        means = pd.read_csv(dir + f"/{config}/means.csv")
        means = means[means['Time'] <= tmax]
        time = means['Time']
        ax.plot(time, means.iloc[:, 1 + comp_index], label=config_labels[c], color = compartment_colors[comp][c])

    ax.set_xlabel("Time [days]")
    ax.set_ylabel(f"Mean {compartment_names[comp]}")
    #ax.set_yscale('log')
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(6,6))
    if(config_labels[0] > ""):
        ax.legend()
    fig.tight_layout()
    fig.savefig(save_dir + f"/all_means_{compartment_names[comp]}_r{region}.png", dpi=dpi)
    plt.close(fig)
    
def plot_mean_errors(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, save_dir, comp_index, figsize):
    fig_mean, ax_mean = plt.subplots(figsize=figsize)
    fig_max, ax_max = plt.subplots(figsize=figsize)
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    errors_mean = []
    errors_max = []
    for c, config in enumerate(configs):        
        means_model1 = pd.read_csv(dir_model1 + f"/{config}/{condition}/{closure_method}/{closure_order}/{init_time}/means.csv")        
        means_model2 = pd.read_csv(dir_model2 + f"/{config}/means.csv")
        err_max = np.max(abs(means_model1.iloc[:, 1 + comp_index] - means_model2.iloc[:, 1 + comp_index])/1.)
        err_mean = np.mean(np.square(abs(means_model1.iloc[:, 1 + comp_index] - means_model2.iloc[:, 1 + comp_index])/1.))
        errors_mean.append(err_mean)
        errors_max.append(err_max)
    
    x = range(len(configs))
    ax_mean.bar(x, errors_mean, color=[colors["dark red"] for _ in range(len(configs))], width=0.5)
    ax_mean.set_xticks(x)
    ax_mean.set_xticklabels(config_labels)
    ax_max.bar(x, errors_max, color=[colors["dark red"] for _ in range(len(configs))], width=0.5)
    ax_max.set_xticks(x)
    ax_max.set_xticklabels(config_labels)
    ax_mean.set_ylabel(r"MSE($\mu_I$)")
    ax_mean.set_yscale("log")
    ax_max.set_ylabel(r"$Max_{t}$ AE($\mu_I$)")
    ax_max.set_yscale("log")
    fig_mean.tight_layout()
    fig_max.tight_layout()
    fig_mean.savefig(save_dir + f"/MSE_mean_{comp}.png", dpi=dpi)
    fig_max.savefig(save_dir + f"/MAX_mean_{comp}.png", dpi=dpi)
    plt.close(fig_mean)
    plt.close(fig_max)
    
def plot_variance_errors(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, save_dir, comp_index, figsize):
    fig_mean, ax_mean = plt.subplots(figsize=figsize)
    fig_max, ax_max = plt.subplots(figsize=figsize)
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    errors_mean = []
    errors_max = []
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        moments_model1 = pd.read_csv(dir_model1 + f"/{config}/{condition}/{closure_method}/{closure_order}/{init_time}/moments.csv")        
        moments_model2 = pd.read_csv(dir_model2 + f"/{config}/moments.csv")
        err_mean = np.mean(np.square(abs(moments_model1[col_name] - moments_model2[col_name])/1.))
        err_max = np.max(abs(moments_model1[col_name] - moments_model2[col_name])/1.)
        errors_mean.append(err_mean)
        errors_max.append(err_max)
    
    x = range(len(configs))
    ax_mean.bar(x, errors_mean, color=[colors["brown"] for _ in range(len(configs))], width=0.5)
    ax_mean.set_xticks(x)
    ax_mean.set_xticklabels(config_labels)
    ax_max.bar(x, errors_max, color=[colors["brown"] for _ in range(len(configs))], width=0.5)
    ax_max.set_xticks(x)
    ax_max.set_xticklabels(config_labels)
    ax_mean.set_ylabel(r"MSE($\sigma^2_I$)")
    ax_mean.set_yscale("log")
    ax_max.set_ylabel(r"$Max_{t}$ AE($\sigma^2_I$)")
    ax_max.set_yscale("log")
    fig_mean.tight_layout()
    fig_max.tight_layout()
    fig_mean.savefig(save_dir + f"/MSE_var_{comp}.png", dpi=dpi)
    fig_max.savefig(save_dir + f"/MAX_var_{comp}.png", dpi=dpi)
    plt.close(fig_mean)
    plt.close(fig_max)
        

def print_max_comp_mean(dir, configs, comp_index, tmax):
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    for c, config in enumerate(configs):        
        means = pd.read_csv(dir + f"/{config}/means.csv")
        means = means[means['Time'] <= tmax]
        means = means.iloc[:, 1 + comp_index]
        print(f"{config}: Max value for mean {compartment_names[comp]} in region {region}: {np.max(means)}")
        
def print_max_comp_variance(dir, configs, comp_index, tmax, num_regions):
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2  # Variance moments: M200, M020, M002
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        moments = pd.read_csv(dir + f"/{config}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        moments = moments[col_name]
        print(f"{config}: Max value for var {compartment_names[comp]} in region {region}: {np.max(moments)}")
        print(f"{config}: Max value for stddev {compartment_names[comp]} in region {region}: {np.max(np.sqrt(moments))}")
        
def print_max_stddev_mean_relation(dir, configs, comp_index, tmax, num_regions):
    region = int(comp_index / len(compartment_names))
    comp = list(compartment_colors.keys())[comp_index % len(compartment_names)]
    moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
    moment_index[comp_index] = 2
    col_name = f"M"
    for m in moment_index:
        col_name += f"{m}"
    for c, config in enumerate(configs):        
        means = pd.read_csv(dir + f"/{config}/means.csv")
        means = means[means['Time'] <= tmax]
        means = means.iloc[:, 1 + comp_index]
        moments = pd.read_csv(dir + f"/{config}/moments.csv")
        moments = moments[moments['Time'] <= tmax]
        moments = moments[col_name]
        print(f"{config}: Max relation value for {compartment_names[comp]} in region {region}: {np.max(np.sqrt(moments)/means)}")

if __name__ == "__main__":
    figsize1 = (4., 2.9)
    figsize2 =(4., 3.2)
    dir = "V:/bick_ju/TemporalHybrid/SMM/SIRS"
    save_dir = "H:/Documents/TemporalHybridModel/SMM/SIRS"
    model2 = "SMM"
    model1 = "Moments"
    model2_name = "Stochastic"
    model1_name = "MoM"
    model2_colors = [colors["very dark blue"], colors["very dark green"], colors["dark teal"], colors["black"], colors["purple"]]
    model1_colors = [colors["middle blue"], colors["middle green"], colors["teal"], colors["middle grey"], colors["rose"]]
    condition = ""
    closure_method = "truncation"
    closure_order = "closure_order_3"
    init_time = "0.000000"
    configs = ["config1_1r_I0_1", "config1_1r_I0_10", "config1_1r_I0_100"]
    config_labels = [r"$I_0=1$", r"$I_0=10$", r"$I_0=100$"]
    num_regions = 1
    tmax = 300
    
    
    
    dir_model1 = f"{dir}/{model1}"
    dir_model2 = f"{dir}/{model2}"
    save_dir = f"{save_dir}"
    os.makedirs(save_dir, exist_ok=True)
    
    # plot_variance_and_std_two_models(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, model1_name, model2_name, save_dir, 1, figsize2, tmax, num_regions, model1_colors, model2_colors)
    # plot_means_two_models(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, model1_name, model2_name, save_dir, 1, figsize2, tmax, num_regions, model1_colors, model2_colors)
    # plot_mean_errors(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, save_dir, 1, figsize1)
    # plot_variance_errors(dir_model1, dir_model2, condition, closure_method, closure_order, init_time, configs, config_labels, save_dir, 1, figsize1)
    
    # print_max_comp_mean(dir, configs, 1, tmax)
    # print_max_comp_variance(dir, configs, 1, tmax, num_regions)
    
    # print_max_stddev_mean_relation(dir, configs, 1, tmax, num_regions)
    
    for r in range(num_regions):
        for comp_index in range(1, 2):
            plot_means(dir, configs, config_labels, save_dir, comp_index + r * len(compartment_names), figsize1, tmax, num_regions)
            plot_variance_and_std(dir, configs, config_labels, save_dir, comp_index + r * len(compartment_names), figsize1, tmax, num_regions)
