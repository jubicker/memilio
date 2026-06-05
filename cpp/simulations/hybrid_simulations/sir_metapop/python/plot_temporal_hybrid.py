from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

def plot_mean_stddev_relation(dir, save_dir, comp_index, num_regions, figsize, tmax):
    comp = list(compartment_colors.keys())[comp_index]
    for r in range(num_regions):
        mean = pd.read_csv(dir + f"/means.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
              
        var = pd.read_csv(dir + f"/moments.csv")
        var = var[var['Time'] <= tmax]
        
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(var.Time, np.sqrt(var[col_name].iloc[:]) / mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color="blue")
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\frac{\sigma_I}{\mu_I}$")
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/relation_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
        

def plot_mean_ts(dir, save_dir, comp_index, num_regions, figsize, tmax, tmin, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            print("Time SMM: ", mean.Time)
            mean.Time.to_csv("smm_time.csv", index=False)
            print("Mean SMM: ", mean.iloc[:, 1 + comp_index + r * len(compartment_names)])
            mean.iloc[:, 1 + comp_index + r * len(compartment_names)].to_csv("smm_mean.csv", index=False)
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=color_smm, label = r"SMM")
        
        # Plot time frame modeled with SMM
        mean = pd.read_csv(dir + f"/means.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        print("Time hybrid SMM: ", mean.Time)
        mean.Time.to_csv("hybrid_smm_time.csv", index=False)
        print("Mean hybrid SMM: ", mean.iloc[:, 1 + comp_index + r * len(compartment_names)])
        mean.iloc[:, 1 + comp_index + r * len(compartment_names)].to_csv("hybrid_smm_mean.csv", index=False)
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=compartment_colors[comp][0], label = r"THMM")
        
        mean = pd.read_csv(dir + f"/expected_values_moments.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        print("Time hybrid moments: ", mean.Time)
        mean.Time.to_csv("hybrid_moments_time.csv", index=False)
        print("Mean hybrid moments: ", mean.iloc[:, 1 + comp_index + r * len(compartment_names)])
        mean.iloc[:, 1 + comp_index + r * len(compartment_names)].to_csv("hybrid_moments_mean.csv", index=False)
        #ax.plot(mean.Time.iloc[1:], mean.iloc[1:, 1 + comp_index + r * len(compartment_names)], color=compartment_colors[comp][0])
        # Add line at swithcing tp
        if(len(mean.Time) > 1):
            times = mean.Time.values
            # Track continuous segments
            segments = []
            start = times[0]

            for i in range(1, len(times)):
                if abs(times[i] - times[i-1]) > 0.11:
                    # End current segment
                    segments.append((start, times[i-1]))
                    start = times[i]  # Start new segment

            # Add last segment
            segments.append((start, times[-1]))

            # Shade continuous regions (grey)
            for start, end in segments:
                ax.axvspan(start, end, color="grey", alpha=0.3)
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\mu_I$")
        if(dir_smm!=""):
            mean_smm = pd.read_csv(dir_smm + f"/means.csv")
            ax.set_ylim(-0.5*np.max(mean_smm.iloc[1:, 1 + comp_index + r * len(compartment_names)]), 2*np.max(mean_smm.iloc[1:, 1 + comp_index]))
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
        
    # One plot for both regions
    fig, ax = plt.subplots(figsize=figsize)
    region_colors_hybrid = [colors["middle blue"], colors["middle green"]]
    region_colors_smm = [colors["dark blue"], colors["dark green"]]
    for r in range(num_regions):
        
        # Plot whole time course of hybrid model
        mean = pd.read_csv(dir + f"/means.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=region_colors_hybrid[r], label = r"Temporal-hybrid model, region " + str(r+1), linewidth=2)
        
        mean = pd.read_csv(dir + f"/expected_values_moments.csv")
        mean = mean[mean['Time'] <= tmax]
        mean = mean[mean['Time'] >= tmin]
        # Add line at swithcing tp
        if(len(mean.Time) > 1):
            times = mean.Time.values
            # Track continuous segments
            segments = []
            start = times[0]

            for i in range(1, len(times)):
                if abs(times[i] - times[i-1]) > 0.11:
                    # End current segment
                    segments.append((start, times[i-1]))
                    start = times[i]  # Start new segment

            # Add last segment
            segments.append((start, times[-1]))

            # Shade continuous regions (grey)
            for start, end in segments:
                ax.axvspan(start, end, color="grey", alpha=0.3)
            
        if(dir_smm) != "":
            mean = pd.read_csv(dir_smm + f"/means.csv")
            mean = mean[mean['Time'] <= tmax]
            mean = mean[mean['Time'] >= tmin]
            ax.plot(mean.Time, mean.iloc[:, 1 + comp_index + r * len(compartment_names)], color=region_colors_smm[r], label = r"Stochastic model, region " + str(r+1), linestyle="dotted", linewidth=2)
            
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(r"$\mu_I$")
    #ax.legend()
    fig.tight_layout()
    fig.savefig(f"{save_dir}/mean_{compartment_names[comp]}_all_regions.png", dpi=dpi)
    plt.close(fig)
        
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/all_regions_legend.png", dpi=dpi, bbox_inches='tight')
        
def plot_var_ts(dir, save_dir, comp_index, num_regions, figsize, tmax, dir_smm="", color_smm=""):
    comp = list(compartment_colors.keys())[comp_index]
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
            
        fig, ax = plt.subplots(figsize=figsize)        
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            var = var[var['Time'] <= tmax]
            var.Time.to_csv("smm_var_time.csv", index=False)
            var[col_name].to_csv("smm_var.csv", index=False)
            ax.plot(var.Time, var[col_name].iloc[:], color=color_smm, label = r"SMM")
        
        # Plot variance
        var = pd.read_csv(dir + f"/moments.csv")
        var = var[var['Time'] <= tmax]
        var.Time.to_csv("hybrid_smm_var_time.csv", index=False)
        var[col_name].to_csv("hybrid_smm_var.csv", index=False)
        ax.plot(var.Time, var[col_name].iloc[:], color = compartment_colors[comp][0], label = r"THMM")
        
        # Moment variance
        var = pd.read_csv(dir + f"/moments_moments.csv")
        var = var[var['Time'] <= tmax]
        var.Time.to_csv("hybrid_moments_var_time.csv", index=False)
        var[col_name].to_csv("hybrid_moments_var.csv", index=False)
        print(f"Moment var at {var.Time.iloc[-2]}: ", var[col_name].iloc[-2])
        #ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = compartment_colors[comp][0])
        # Add line at swithcing tp
        if(len(var.Time) > 1):
            times = var.Time.values
            # Track continuous segments
            segments = []
            start = times[0]

            for i in range(1, len(times)):
                if abs(times[i] - times[i-1]) > 0.11:
                    # End current segment
                    segments.append((start, times[i-1]))
                    start = times[i]  # Start new segment

            # Add last segment
            segments.append((start, times[-1]))

            # Shade continuous regions (grey)
            for start, end in segments:
                ax.axvspan(start, end, color="grey", alpha=0.3)
            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\sigma^2_I$")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
        #ax.set_yscale("log")
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_r{r}.png", dpi=dpi)
        plt.close(fig)
    
    # All regions in one plot
    region_colors_hybrid = [colors["middle blue"], colors["middle green"]]
    region_colors_smm = [colors["dark blue"], colors["dark green"]]
    fig, ax = plt.subplots(figsize=figsize)    
    for r in range(num_regions):
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        
        # Plot variance
        var = pd.read_csv(dir + f"/moments.csv")
        var = var[var['Time'] <= tmax]
        ax.plot(var.Time, var[col_name].iloc[:], color = region_colors_hybrid[r], label = r"THMM, region " + str(r+1), linewidth=2)
        
        # Moment variance
        var = pd.read_csv(dir + f"/moments_moments.csv")
        var = var[var['Time'] <= tmax]
        #ax.plot(var.Time.iloc[1:], var[col_name].iloc[1:], color = region_colors_hybrid[r], linewidth=2)
        # Add line at swithcing tp
        if(len(var.Time) > 1):
            times = var.Time.values
            # Track continuous segments
            segments = []
            start = times[0]

            for i in range(1, len(times)):
                if abs(times[i] - times[i-1]) > 0.11:
                    # End current segment
                    segments.append((start, times[i-1]))
                    start = times[i]  # Start new segment

            # Add last segment
            segments.append((start, times[-1]))

            # Shade continuous regions (grey)
            for start, end in segments:
                ax.axvspan(start, end, color="grey", alpha=0.3)

            
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"$\sigma^2_I$")
        if(dir_smm!=""):
            var_smm = pd.read_csv(dir_smm + f"/moments.csv")
            ax.set_ylim(-0.5*np.max(var_smm[col_name].iloc[:]), 1.5*np.max(var_smm[col_name].iloc[:]))
            
        if(dir_smm) != "":
            var = pd.read_csv(dir_smm + f"/moments.csv")
            var = var[var['Time'] <= tmax]
            ax.plot(var.Time, var[col_name].iloc[:], color=region_colors_smm[r], label = r"SMM, region " + str(r+1), linestyle="dotted", linewidth=2)
        #ax.set_yscale("log")
        #ax.legend()
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_{compartment_names[comp]}_all_regions.png", dpi=dpi)
        plt.close(fig)

def plot_mean_error(dir, dir_smm, figsize, save_dir, num_regions, comp_index, tmax):
    fig, ax = plt.subplots(figsize=figsize)
    errors_mean = []
    errors_max = []
    for r in range(num_regions):
        mean_smm = pd.read_csv(dir_smm + f"/means.csv")
        mean_smm = mean_smm[mean_smm.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)]
        mean_hybrid = pd.read_csv(dir + f"/means.csv")
        mean_hybrid = mean_hybrid[mean_hybrid.Time <= tmax].iloc[:, 1 + comp_index + r * len(compartment_names)]
        err =  np.mean(np.square(abs(mean_hybrid - mean_smm)/1.))
        errors_mean.append(err)
        err_max = np.max(abs(mean_hybrid - mean_smm)/1.)
        errors_max.append(err_max)
    x = np.arange(num_regions)
    bar_width = 0.4
    ax.bar(x - bar_width/2, errors_mean, color=colors['dark red'], width=bar_width, label=r"MSE")
    ax.bar(x + bar_width/2, errors_max, color=colors['brown'], width=bar_width, label=r"Max AE")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Region " + str(r+1) for r in range(num_regions)])
    ax.set_ylabel(r"err($\mu_I$)")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(f"{save_dir}/error_all_regions.png", dpi=dpi)
    plt.close(fig)
    
    fig_leg = plt.figure(figsize=(6, 2))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    handles, labels = ax.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc='center', ncol=1)

    fig_leg.savefig(save_dir + "/error_legend.png", dpi=dpi, bbox_inches='tight')

def plot_switch_tp_vs_error_mean(dir, dir_smm, ode_dir, conditions, closure_order, closure_method, figsize, save_dir, num_regions, comp_index, condition_names):
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        fig_ts, ax_ts = plt.subplots(figsize=figsize)
        ax2 = ax.twinx()
        errors_mean = []
        errors_max = []
        switch_tps = []
        mean_smm = pd.read_csv(dir_smm + f"/means.csv").iloc[:, 1 + comp_index + r * len(compartment_names)]
        mean_ode = pd.read_csv(ode_dir + f"/{closure_method}/{closure_order}/0.000000/means.csv").iloc[:, 1 + comp_index + r * len(compartment_names)]
        err_ode = np.mean(np.square(abs(mean_ode - mean_smm)/1.))
        ax_ts.plot(abs(mean_ode - mean_smm), color=colors['dark red'], label="ODE")
        errors_mean.append(err_ode)
        errors_max.append(np.max(abs(mean_ode - mean_smm)/1.))
        switch_tps.append(0.)
        for i, condition in enumerate(conditions):
            mean_hybrid = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/means.csv").iloc[:, 1 + comp_index + r * len(compartment_names)]
            err =  np.mean(np.square(abs(mean_hybrid - mean_smm)/1.))
            ax_ts.plot(abs(mean_hybrid - mean_smm), label=condition, color=list(colors.values())[i])
            errors_mean.append(err)
            err_max = np.max(abs(mean_hybrid - mean_smm)/1.)
            errors_max.append(err_max)
            switch_tp = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/expected_values_moments.csv").Time.iloc[1]
            switch_tps.append(switch_tp)
        x = np.arange(len(condition_names) + 1)
        bar_width = 0.4
        ax.bar(x - bar_width/2, errors_mean, color=colors['dark red'], width=bar_width)
        ax.bar(x + bar_width/2, errors_max, color=colors['brown'], width=bar_width)
        ax2.plot(x, switch_tps, color=colors['middle blue'], marker='o')
        ax.set_xticks(x)
        ax.set_xticklabels(["MoM"] + condition_names)
        ax.set_ylabel(r"err($\mu_I$)")
        ax2.set_ylabel("Switching time [days]")
        ax.set_yscale("log")
        fig.tight_layout()
        fig.savefig(f"{save_dir}/error_switch_tp_r{r}.png", dpi=dpi)
        ax_ts.set_yscale("log")
        fig_ts.savefig(f"{save_dir}/error_ts_r{r}.png", dpi=dpi)
        plt.close(fig)
        
def plot_switch_tp_vs_error_var(dir, dir_smm, ode_dir, conditions, closure_order, closure_method, figsize, save_dir, num_regions, comp_index, condition_names):
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        fig_ts, ax_ts = plt.subplots(figsize=figsize)
        ax2 = ax.twinx()
        errors_mean = []
        errors_max = []
        switch_tps = []
        moment_index = [0 for _ in range(len(compartment_names)*num_regions)]
        moment_index[comp_index + r * len(compartment_names)] = 2
        col_name = f"M"
        for m in moment_index:
            col_name += f"{m}"
        mean_smm = pd.read_csv(dir_smm + f"/moments.csv")[col_name].iloc[:]
        mean_ode = pd.read_csv(ode_dir + f"/{closure_method}/{closure_order}/0.000000/moments.csv")[col_name].iloc[:]
        err_ode = np.mean(np.square(abs(mean_ode - mean_smm)/1.))
        ax_ts.plot(abs(mean_ode - mean_smm), color=colors['dark red'], label="ODE")
        errors_mean.append(err_ode)
        errors_max.append(np.max(abs(mean_ode - mean_smm)/1.))
        switch_tps.append(0.)
        for i, condition in enumerate(conditions):
            mean_hybrid = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/moments.csv")[col_name].iloc[:]
            err =  np.mean(np.square(abs(mean_hybrid - mean_smm)/1.))
            ax_ts.plot(abs(mean_hybrid - mean_smm), label=condition, color=list(colors.values())[i])
            errors_mean.append(err)
            err_max = np.max(abs(mean_hybrid - mean_smm)/1.)
            errors_max.append(err_max)
            switch_tp = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/expected_values_moments.csv").Time.iloc[1]
            switch_tps.append(switch_tp)
        x = np.arange(len(condition_names) + 1)
        bar_width = 0.4
        ax.bar(x - bar_width/2, errors_mean, color=colors['dark red'], width=bar_width)
        ax.bar(x + bar_width/2, errors_max, color=colors['brown'], width=bar_width)
        ax2.plot(x, switch_tps, color=colors['middle blue'], marker='o')
        ax.set_xticks(x)
        ax.set_xticklabels(["MoM"] + condition_names)
        ax.set_ylabel(r"err($\sigma_I^2$)")
        ax.set_yscale("log")
        ax2.set_ylabel("Switching time [days]")
        fig.tight_layout()
        fig.savefig(f"{save_dir}/var_error_switch_tp_r{r}.png", dpi=dpi)
        ax_ts.set_yscale("log")
        fig_ts.savefig(f"{save_dir}/var_error_ts_r{r}.png", dpi=dpi)
        plt.close(fig)
        
def plot_switch_tp_vs_runtime(dir, dir_smm, ode_dir, conditions, closure_order, closure_method, figsize, save_dir, num_regions, comp_index, condition_names):
    for r in range(num_regions):
        fig, ax = plt.subplots(figsize=figsize)
        ax2 = ax.twinx()
        runtimes = []
        switch_tps = []
        runtime_smm = pd.read_csv(dir_smm + f"/total_time.csv").Runtime.iloc[0]
        runtime_ode = pd.read_csv(ode_dir + f"/{closure_method}/{closure_order}/0.000000/total_time.csv").Runtime.iloc[0]
        runtimes.append(runtime_ode)
        switch_tps.append(0.)
        for i, condition in enumerate(conditions):
            runtime_hybrid = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/total_time.csv").Runtime.iloc[0]
            runtimes.append(runtime_hybrid)
            switch_tp = pd.read_csv(dir + f"/{condition}/{closure_method}/{closure_order}/expected_values_moments.csv").Time.iloc[1]
            switch_tps.append(switch_tp)
        runtimes.append(runtime_smm)
        x = range(len(condition_names) + 2)
        ax.bar(x, runtimes, color=colors['purple'], width=0.5)
        ax.set_xticks(x)
        x = range(len(condition_names) + 1)
        ax2.plot(x, switch_tps, color=colors['middle blue'], marker='o')
        ax.set_xticklabels(["MoM"] + condition_names + ["Stoch"])
        ax.set_ylabel(r"Runtime [s]")
        ax.set_yscale("log")
        ax2.set_ylabel("Switching time [days]")
        fig.tight_layout()
        fig.savefig(f"{save_dir}/runtime_switch_tp_r{r}.png", dpi=dpi)
        plt.close(fig)
    
    #Legend    
    legend_elements = [
    Patch(facecolor=colors['purple'], edgecolor=colors['purple'], label='Runtime'),
    Patch(facecolor=colors['dark red'], edgecolor=colors['dark red'], label='MSE'),
    Patch(facecolor=colors['brown'], edgecolor=colors['brown'], label=r'$Max_t$ AE'),
    Line2D([0], [0], color=colors['middle blue'], marker='o', label='Switching time')
    ]
    fig, ax = plt.subplots()
    ax.legend(handles=legend_elements, loc='center')
    ax.axis('off')
    plt.savefig(f"{save_dir}/legend.png", bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    figsize = (4, 3)
    dir = "V:/bick_ju/TemporalHybrid"
    save_dir = "H:/Documents/TemporalHybridModel"
    hybrid_model = "Hybrid2"
    config = "config1_1r_I0_1000"
    num_regions = 1
    condition = "mean_stddev_relation"
    conditions = ["mean_threshold/0.001000", "mean_threshold/0.010000/", "mean_threshold/0.100000", "mean_threshold/0.300000", "var_gradient"]
    condition_names = [r"$\tau_{\mu_I}=$" +"\n" + r"$0.1$%", r"$\tau_{\mu_I}=$" +"\n" + r"$1$%", r"$\tau_{\mu_I}=$" + "\n" + r"$10$%", r"$\tau_{\mu_I}=$" + "\n" + r"$30$%", r"$\nabla \sigma_I^2$"]
    closure_method = "truncation"
    closure_order = "closure_order_3"
    color_smm = colors['dark grey']
    tmin = 0
    tmax = 90
    condition_name = ""
    
    smm_dir = f"{dir}/SMM/{config}"
    hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
    os.makedirs(save_dir, exist_ok=True)
    
    comp_index = 1
    plot_mean_ts(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax, tmin, smm_dir, color_smm)
    plot_var_ts(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax=tmax, dir_smm=smm_dir, color_smm=color_smm)
    plot_mean_stddev_relation(hybrid_dir, save_dir, comp_index, num_regions, figsize, tmax)
    #plot_mean_error(hybrid_dir, smm_dir, figsize, save_dir, num_regions, comp_index, tmax)
    
    # plot_switch_tp_vs_error_mean(f"{dir}/{hybrid_model}/{config}", smm_dir, ode_dir=f"{dir}/Moments/{config}", conditions=conditions, closure_order=closure_order, closure_method=closure_method, figsize=figsize, save_dir="H:/Documents/TemporalHybridModel" + f"/{hybrid_model}/{config}", num_regions=num_regions, comp_index=comp_index, condition_names=condition_names)
    # plot_switch_tp_vs_error_var(f"{dir}/{hybrid_model}/{config}", smm_dir, ode_dir=f"{dir}/Moments/{config}", conditions=conditions, closure_order=closure_order, closure_method=closure_method, figsize=figsize, save_dir="H:/Documents/TemporalHybridModel" + f"/{hybrid_model}/{config}", num_regions=num_regions, comp_index=comp_index, condition_names=condition_names)
    # plot_switch_tp_vs_runtime(f"{dir}/{hybrid_model}/{config}", smm_dir, ode_dir=f"{dir}/Moments/{config}", conditions=conditions, closure_order=closure_order, closure_method=closure_method, figsize=figsize, save_dir="H:/Documents/TemporalHybridModel" + f"/{hybrid_model}/{config}", num_regions=num_regions, comp_index=comp_index, condition_names=condition_names)
