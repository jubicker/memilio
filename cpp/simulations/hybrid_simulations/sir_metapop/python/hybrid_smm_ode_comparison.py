from settings import *
import os

def plot_means(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/means.csv")
        figsize = (6,4)
        time = means["Time"]
        ax_S.plot(time, means['muS_r0'], label = "Switch" + sv)
        ax_I.plot(time, means['muI_r0'], label = "Switch" + sv)
        ax_R.plot(time, means['muR_r0'], label = "Switch" + sv)
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Susceptible [#]")
    ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "mean_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Infected [#]")
    ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "mean_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Recovered [#]")
    ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "mean_R.png", dpi=dpi)
    plt.close(fig_R)
    
def plot_variances(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        figsize = (6,4)
        time = means["Time"]
        ax_S.plot(time, means['M200'], label = "Switch" + sv)
        ax_I.plot(time, means['M020'], label = "Switch" + sv)
        ax_R.plot(time, means['M002'], label = "Switch" + sv)
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Var(S)")
    ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "var_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Var(I)")
    ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "var_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Var(R)")
    ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "var_R.png", dpi=dpi)
    plt.close(fig_R)

def plot_covariances(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_1, ax_1 = plt.subplots(figsize=figsize)
    fig_2, ax_2 = plt.subplots(figsize=figsize)
    fig_3, ax_3 = plt.subplots(figsize=figsize)
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        figsize = (6,4)
        time = means["Time"]
        ax_1.plot(time, means['M110'], label = "Switch" + sv)
        ax_2.plot(time, means['M101'], label = "Switch" + sv)
        ax_3.plot(time, means['M011'], label = "Switch" + sv)
        
    ax_1.set_xlabel("Time [days]")
    ax_1.set_ylabel("Cov(S,I)")
    ax_1.legend()
    fig_1.tight_layout()
    fig_1.savefig(save_dir + "cov_SI.png", dpi=dpi)
    plt.close(fig_1)
    
    ax_2.set_xlabel("Time [days]")
    ax_2.set_ylabel("Cov(S,R)")
    ax_2.legend()
    fig_2.tight_layout()
    fig_2.savefig(save_dir + "cov_SR.png", dpi=dpi)
    plt.close(fig_2)
    
    ax_3.set_xlabel("Time [days]")
    ax_3.set_ylabel("Cov(R,I)")
    ax_3.legend()
    fig_3.tight_layout()
    fig_3.savefig(save_dir + "cov_RI.png", dpi=dpi)
    plt.close(fig_3)

def plot_err_mean(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    
    MAEs = {}
    RMSEs = {}
    MAPEs = {}
    xlabels = []
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/means.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
    
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0] 
        ax_S.plot(time, (df['muS_r0'] - true['muS_r0']).abs(), label = "Switch" + sv)
        ax_I.plot(time, (df['muI_r0'] - true['muI_r0']).abs(), label = "Switch" + sv)
        ax_R.plot(time, (df['muR_r0'] - true['muR_r0']).abs(), label = "Switch" + sv)
        xlabels.append(sv[:6])
        
        # Calculate total error
        MAEs[sv] = {
            'S': (df['muS_r0'] - true['muS_r0']).abs().mean(),
            'I': (df['muI_r0'] - true['muI_r0']).abs().mean(),
            'R': (df['muR_r0'] - true['muR_r0']).abs().mean()
        }
        RMSEs[sv] = {
            'S': np.sqrt(((df['muS_r0'] - true['muS_r0'])**2).mean()),
            'I': np.sqrt(((df['muI_r0'] - true['muI_r0'])**2).mean()),
            'R': np.sqrt(((df['muR_r0'] - true['muR_r0'])**2).mean())
        }
        MAPEs[sv] = {
            'S': ( (df['muS_r0'] - true['muS_r0']).abs() / true['muS_r0'].replace(0, np.nan) ).mean(),
            'I': ( (df['muI_r0'] - true['muI_r0']).abs() / true['muI_r0'].replace(0, np.nan) ).mean(),
            'R': ( (df['muR_r0'] - true['muR_r0']).abs() / true['muR_r0'].replace(0, np.nan) ).mean()
        }
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Err(S)")
    ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Err(I)")
    ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Err(R)")
    ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_R.png", dpi=dpi)
    plt.close(fig_R)
    
        # Grouped bar plot with total MAE and MSE
    labels = list(MAEs.keys())
    x = range(len(labels))
    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=(6, 4))
    ax_bar_S_MAPE = ax_bar_S.twinx()
    vals_mae_S = [MAEs[k]['S'] for k in labels]
    vals_rmse_S = [RMSEs[k]['S'] for k in labels]
    h_mae_S = ax_bar_S.bar([i - 0.2 for i in x], vals_mae_S, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_S = ax_bar_S.bar([i + 0.2 for i in x], vals_rmse_S, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_S = ax_bar_S_MAPE.bar(x, [MAPEs[k]['S'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(labels)
    ax_bar_S.set_yscale('log')
    ax_bar_S_MAPE.set_yscale('log')
    ax_bar_S_MAPE.set_ylabel("MAPE")
    ax_bar_S.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_S.legend([h_mae_S[0], h_rmse_S[0], h_mape_S[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Err_S.png", dpi=dpi)
    plt.close(fig_bar_S)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=(6, 4))
    ax_bar_I_MAPE = ax_bar_I.twinx()
    vals_mae_I = [MAEs[k]['I'] for k in labels]
    vals_rmse_I = [RMSEs[k]['I'] for k in labels]
    h_mae_I = ax_bar_I.bar([i - 0.2 for i in x], vals_mae_I, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_I = ax_bar_I.bar([i + 0.2 for i in x], vals_rmse_I, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_I = ax_bar_I_MAPE.bar(x, [MAPEs[k]['I'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(labels)
    ax_bar_I.set_yscale('log')
    ax_bar_I_MAPE.set_yscale('log')
    ax_bar_I_MAPE.set_ylabel("MAPE")
    ax_bar_I.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_I.legend([h_mae_I[0], h_rmse_I[0], h_mape_I[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Err_I.png", dpi=dpi)
    plt.close(fig_bar_I)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=(6, 4))
    ax_bar_R_MAPE = ax_bar_R.twinx()
    vals_mae_R = [MAEs[k]['R'] for k in labels]
    vals_rmse_R = [RMSEs[k]['R'] for k in labels]
    h_mae_R = ax_bar_R.bar([i - 0.2 for i in x], vals_mae_R, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_R = ax_bar_R.bar([i + 0.2 for i in x], vals_rmse_R, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_R = ax_bar_R_MAPE.bar(x, [MAPEs[k]['R'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(labels)
    ax_bar_R.set_yscale('log')
    ax_bar_R_MAPE.set_yscale('log')
    ax_bar_R_MAPE.set_ylabel("MAPE")
    ax_bar_R.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_R.legend([h_mae_R[0], h_rmse_R[0], h_mape_R[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Err_R.png", dpi=dpi)
    plt.close(fig_bar_R)
        
def plot_err_var(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    
    MAEs = {}
    RMSEs = {}
    MAPEs = {}  # added mapes
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
    
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0] 
        ax_S.plot(time, (df['M200'] - true['M200']).abs(), label = "Switch" + sv)
        ax_I.plot(time, (df['M020'] - true['M020']).abs(), label = "Switch" + sv)
        ax_R.plot(time, (df['M002'] - true['M002']).abs(), label = "Switch" + sv)
        
        # Calculate total error
        MAEs[sv] = {
            'S': (df['M200'] - true['M200']).abs().mean(),
            'I': (df['M020'] - true['M020']).abs().mean(),
            'R': (df['M002'] - true['M002']).abs().mean()
        }
        RMSEs[sv] = {
            'S': np.sqrt(((df['M200'] - true['M200'])**2).mean()),
            'I': np.sqrt(((df['M020'] - true['M020'])**2).mean()),
            'R': np.sqrt(((df['M002'] - true['M002'])**2).mean())
        }
        MAPEs[sv] = {
            'S': ( (df['M200'] - true['M200']).abs() / true['M200'].replace(0, np.nan) ).mean(),
            'I': ( (df['M020'] - true['M020']).abs() / true['M020'].replace(0, np.nan) ).mean(),
            'R': ( (df['M002'] - true['M002']).abs() / true['M002'].replace(0, np.nan) ).mean()
        }
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Err(Var(S))")
    ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_VarS.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Err(Var(I))")
    ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_VarI.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Err(Var(R))")
    ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_VarR.png", dpi=dpi)
    plt.close(fig_R)
    
    # Grouped bar plot with total MAE, RMSE and MAPE
    labels = list(MAEs.keys())
    x = range(len(labels))
    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=(6, 4))
    ax_bar_S_MAPE = ax_bar_S.twinx()
    vals_mae_S = [MAEs[k]['S'] for k in labels]
    vals_rmse_S = [RMSEs[k]['S'] for k in labels]
    h_mae_S = ax_bar_S.bar([i - 0.2 for i in x], vals_mae_S, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_S = ax_bar_S.bar([i + 0.2 for i in x], vals_rmse_S, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_S = ax_bar_S_MAPE.bar(x, [MAPEs[k]['S'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(labels)
    ax_bar_S.set_yscale('log')
    ax_bar_S_MAPE.set_yscale('log')
    ax_bar_S_MAPE.set_ylabel("MAPE")
    ax_bar_S.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_S.legend([h_mae_S[0], h_rmse_S[0], h_mape_S[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Err_VarS.png", dpi=dpi)
    plt.close(fig_bar_S)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=(6, 4))
    ax_bar_I_MAPE = ax_bar_I.twinx()
    vals_mae_I = [MAEs[k]['I'] for k in labels]
    vals_rmse_I = [RMSEs[k]['I'] for k in labels]
    h_mae_I = ax_bar_I.bar([i - 0.2 for i in x], vals_mae_I, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_I = ax_bar_I.bar([i + 0.2 for i in x], vals_rmse_I, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_I = ax_bar_I_MAPE.bar(x, [MAPEs[k]['I'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(labels)
    ax_bar_I.set_yscale('log')
    ax_bar_I_MAPE.set_yscale('log')
    ax_bar_I_MAPE.set_ylabel("MAPE")
    ax_bar_I.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_I.legend([h_mae_I[0], h_rmse_I[0], h_mape_I[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Err_VarI.png", dpi=dpi)
    plt.close(fig_bar_I)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=(6, 4))
    ax_bar_R_MAPE = ax_bar_R.twinx()
    vals_mae_R = [MAEs[k]['R'] for k in labels]
    vals_rmse_R = [RMSEs[k]['R'] for k in labels]
    h_mae_R = ax_bar_R.bar([i - 0.2 for i in x], vals_mae_R, width=0.2, label='MAE', color = colors['dark blue'])
    h_rmse_R = ax_bar_R.bar([i + 0.2 for i in x], vals_rmse_R, width=0.2, label='RMSE', color = colors['dark green'])
    h_mape_R = ax_bar_R_MAPE.bar(x, [MAPEs[k]['R'] for k in labels], width=0.2, label='MAPE', color = colors['brown'])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(labels)
    ax_bar_R.set_yscale('log')
    ax_bar_R_MAPE.set_yscale('log')
    ax_bar_R_MAPE.set_ylabel("MAPE")
    ax_bar_R.set_ylabel("MAE / RMSE")
    # combined legend
    fig_bar_R.legend([h_mae_R[0], h_rmse_R[0], h_mape_R[0]], ['MAE','RMSE','MAPE'], bbox_to_anchor = (0.8, 0.9))
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Err_VarR.png", dpi=dpi)
    plt.close(fig_bar_R)
    
def plot_runtimes(result_dir, save_dir, switching_values):
    mean_runtimes = []
    median_runtimes = []
    xlabels = []
    for sv in switching_values:
        runtime = pd.read_csv(result_dir + f"switch_value_{sv}/runtimes.csv")
        mean_runtimes.append(runtime['runtime'].mean())
        median_runtimes.append(runtime['runtime'].median())
        xlabels.append(sv[:6])
    
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(xlabels, mean_runtimes, color=colors['dark red'])
    ax.set_xlabel("Switching Value")
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale('log')
    fig.tight_layout()
    fig.savefig(save_dir + "runtimes_mean.png", dpi=dpi)
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(xlabels, median_runtimes, color=colors['dark red'])
    ax.set_xlabel("Switching Value")
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale('log')
    fig.tight_layout()
    fig.savefig(save_dir + "runtimes_median.png", dpi=dpi)
    plt.close(fig)

config = "config_1r"         
result_dir = "V:/bick_ju/TemporalHybrid/Hybrid1/" + config + "/"
switching_values = ["0.000000", "0.000100", "0.001000", "0.010000", "0.100000", "0.300000", "1.000000"]
save_dir = "H:/Documents/TemporalHybridModel/Hybrid1/" + config + "/"
os.makedirs(save_dir, exist_ok=True)

# plot_means(result_dir, save_dir, switching_values)
# plot_variances(result_dir, save_dir, switching_values)
# plot_covariances(result_dir, save_dir, switching_values)
plot_err_mean(result_dir, save_dir, switching_values)
plot_err_var(result_dir, save_dir, switching_values)
# plot_runtimes(result_dir, save_dir, switching_values)
