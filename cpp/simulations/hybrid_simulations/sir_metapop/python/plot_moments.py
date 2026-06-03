from settings import *
import os

def plot_mean_ts(dir, save_dir, closures, comp_index, num_regions, figsize, closure_colors, closure_order, init_tp, tmax, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        # Plot smm mean
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            mean = mean[mean['Time'] <= tmax]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_smm, label = "SMM")
        for c in closures:            
            mean = pd.read_csv(dir + f"/{c}/{init_tp}/" + "means.csv")
            mean = mean[mean['Time'] <= tmax]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=closure_colors[c], label = f"{c}")
                
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(f"{compartment_names[comp]} [#]")
        if(dir_smm!=""):
            mean_smm = pd.read_csv(dir_smm + f"/means.csv")
            ax.set_ylim(-1.5*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)]) - 10, 2*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)])+10)
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        
        handles, labels = ax.get_legend_handles_labels()
        fig_leg = plt.figure(figsize=figsize)                   
        fig_leg.legend(handles, labels, loc='center')      
        fig_leg.tight_layout()
        fig_leg.savefig(save_dir + f"/legend_{comp}.png", dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig_leg)
        plt.close(fig)
        
def plot_var_ts(dir, save_dir, closures, comp_index, num_regions, figsize, closure_colors, closure_order, init_tp, tmax, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2  # Variance moments: M200, M020, M002
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        fig, ax = plt.subplots(figsize=figsize)
        # Plot smm mean
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            var = var[var['Time'] <= tmax]
            ax.plot(var.Time, var[col_name], color=color_smm, label = "SMM")
        for c in closures:           
            var = pd.read_csv(dir + f"/{c}/{init_tp}/" + "moments.csv")
            var = var[var['Time'] <= tmax]
            ax.plot(var.Time, var[col_name], color=closure_colors[c], label = f"{c}")
                
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(f"{compartment_names[comp]} [#]")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
        
def plot_mean_error(smm_dir, hybrid_dir, save_dir, closure_names, colors_hybrid, comp_index, num_regions, figsize, closure_order = "closure_order_3"):
    comp = list(compartment_colors.keys())[comp_index]
    fig, ax = plt.subplots(figsize=figsize)
    
    smm_values = pd.read_csv(smm_dir + "/means.csv")
    
    summed_errors = {r: [] for r in range(num_regions)}
    linewidth = 1
    
    for r in range(num_regions):
    
        for i, st in enumerate(closure_names):
            
            values = pd.read_csv(hybrid_dir + f"/{st}/{closure_order}/means.csv")
            error = (smm_values.iloc[:-1, 1 + comp_index + r * len(compartment_names)] - values.iloc[:, 1 + comp_index + r * len(compartment_names)]).abs()
            summed_errors[r].append(error.sum())
            label = f"{st}"
            if(num_regions > 1):
                label += f" Region {r}" 
            ax.plot(values.Time, error, label = label, color = colors_hybrid[r][i], linewidth=linewidth)
        
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"|$\mu_t - \hat{\mu}_t$|")
        ax.set_yscale("log")
        fig.tight_layout()
        fig.savefig(save_dir + f"/mean_error_ts_{comp}_r{r}", dpi=dpi)
    
    # Save legend
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center')      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + f"/legend_{comp}.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    # Save summed error
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(closure_names))
    width = 0.8 / num_regions
    for r in range(num_regions):
        summed_errors_region = summed_errors[r]
        bar_colors = colors_hybrid[r]
        ax.bar([i + (r-1) * width/2 + r * width/2 for i in x], summed_errors_region, width=width, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel(r"$\sum_{t}$|$\mu_t - \hat{\mu}_t$|")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/mean_sum_err_{comp}.png", dpi=dpi)
    plt.close()
        
def plot_var_error(smm_dir, hybrid_dir, save_dir, closure_names, colors_hybrid, comp_index, num_regions, figsize, closure_order = "closure_order_3"):
    comp = list(compartment_colors.keys())[comp_index]
    fig, ax = plt.subplots(figsize=figsize)
    
    smm_values = pd.read_csv(smm_dir + "/moments.csv")
    
    summed_errors = {r: [] for r in range(num_regions)}
    linewidth = 1
    
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2  # Variance moments: M200, M020, M002
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        for i, st in enumerate(closure_names):
            
            values = pd.read_csv(hybrid_dir + f"/{st}/{closure_order}/moments.csv")
            error = (smm_values[col_name].iloc[:-1] - values[col_name].iloc[:]).abs()
            summed_errors[r].append(error.sum())
            label = f"{st}"
            if(num_regions > 1):
                label += f" Region {r}" 
            ax.plot(values.Time, error, label = label, color = colors_hybrid[r][i], linewidth=linewidth)
        
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"|$\sigma^2_t - \hat{\sigma^2}_t$|")
        ax.set_yscale("log")
        fig.tight_layout()
        fig.savefig(save_dir + f"/var_error_ts_{comp}_r{r}", dpi=dpi)
    
    # Save legend
    handles, labels = ax.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center')      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + f"/legend_{comp}.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    # Save summed error
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]*0.8))
    x = range(len(closure_names))
    width = 0.8 / num_regions
    for r in range(num_regions):
        summed_errors_region = summed_errors[r]
        bar_colors = colors_hybrid[r]
        ax.bar([i + (r-1) * width/2 + r * width/2 for i in x], summed_errors_region, width=width, color =  bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel(r"$\sum_{t}$|$\sigma^2_t - \hat{\sigma^2}_t$|")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + f"/var_sum_err_{comp}.png", dpi=dpi)
    plt.close()

if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    config = "config_SIR_I0_0_1_10_100_exchange"
    closure_order = 3
    num_regions = 4
    closures = ["truncation/closure_order_3"]
    start_tp = "0.000000"
    colors_closures = {"truncation": colors["dark blue"], "pair_approx": colors["teal"], "lognorm": colors["dark green"], "lognorm/closure_order_3": colors["middle blue"], "lognorm_zero_infl/closure_order_4": colors["red"], "truncation/closure_order_3": colors["dark green"], "pair_approx/closure_order_3": colors["orange"]}
    color_ode = [colors['black'], colors['dark grey']]
    color_smm = colors['dark grey']
    condition_name = "mean_threshold/0.000000"
    hybrid_model = "Hybrid2"
    compare_values = ["truncation", "pair_approx", "lognorm"]
    colors_hybrid = [[colors['dark blue'], colors['middle blue'], colors['light blue'], colors['teal'], colors['light teal'], colors['dark green'], colors['middle green'], colors['light green']], [colors['purple'], colors['rose'], colors['red'], colors['dark red'], colors['brown']]]
    
    smm_dir = f"{dir}/SMM/{config}"
    moment_dir = f"{dir}/Moments/{config}"
    save_dir = f"{save_dir}/Moments/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition_name}"
    os.makedirs(save_dir, exist_ok=True)
    
    tmax = 3000
    
    comp_index = 1
    
    plot_mean_ts(moment_dir, save_dir, closures, comp_index, num_regions, figsize, colors_closures, closure_order, start_tp, tmax, smm_dir, color_smm)
    plot_var_ts(moment_dir, save_dir, closures, comp_index, num_regions, figsize, colors_closures, closure_order, start_tp, tmax, smm_dir, color_smm)
    
    # plot_mean_error(smm_dir, hybrid_dir, save_dir, compare_values, colors_hybrid, comp_index, num_regions, figsize)
    # plot_var_error(smm_dir, hybrid_dir, save_dir, compare_values, colors_hybrid, comp_index, num_regions, figsize)
