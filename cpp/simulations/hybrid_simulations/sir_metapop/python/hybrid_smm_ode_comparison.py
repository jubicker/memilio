from settings import *
import os

def plot_means(result_dir, save_dir, switching_values, num_regions, my_colors):
    figsize = (3.3, 2.5)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    c = 0
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/means.csv")
        time = means["Time"]
        label = f"{float(sv) * 100:.4f}"
        for r in range(num_regions):            
            ax_S.plot(time, means[f'muS_r{r}'], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c])
            ax_I.plot(time, means[f'muI_r{r}'], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c])
            ax_R.plot(time, means[f'muR_r{r}'], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c])
            c += 1
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Mean(S)")
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "mean_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Mean(I)")
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "mean_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Mean(R)")
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "mean_R.png", dpi=dpi)
    plt.close(fig_R)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
def plot_variances(result_dir, save_dir, switching_values, num_regions, my_colors):
    figsize = (3.3, 2.5)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    c = 0
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        time = means["Time"]
        label = f"{float(sv) * 100:.4f}"
        for r in range(num_regions):
            string_S = "M"
            string_I = "M"
            string_R = "M"
            for _ in range(3*num_regions):
                string_S += "0"
                string_I += "0"
                string_R += "0"
            string_S = string_S[:(r * 3 + 1)] + "2" + string_S[(r * 3 + 1 + 1):]
            string_I = string_I[:(r * 3 + 2)] + "2" + string_I[(r * 3 + 2 + 1):]
            string_R = string_R[:(r * 3 + 3)] + "2" + string_R[(r * 3 + 3 + 1):]
            ax_S.plot(time, means[string_S], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c], linewidth = 5)
            ax_I.plot(time, means[string_I], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c])
            ax_R.plot(time, means[string_R], label = str(label) + f"%, Region {r:.0f}", color = my_colors[c])
            c += 1
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Var(S)")
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "var_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Var(I)")
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "var_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Var(R)")
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "var_R.png", dpi=dpi)
    plt.close(fig_R)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

def plot_covariances(result_dir, save_dir, switching_values):
    figsize = (6,4)
    
    fig_1, ax_1 = plt.subplots(figsize=figsize)
    fig_2, ax_2 = plt.subplots(figsize=figsize)
    fig_3, ax_3 = plt.subplots(figsize=figsize)
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        figsize = (6,4)
        time = means["Time"]
        label = float(sv) * 100.
        ax_1.plot(time, means['M110'], label = "Switch" + str(label) + "%")
        ax_2.plot(time, means['M101'], label = "Switch" + str(label) + "%")
        ax_3.plot(time, means['M011'], label = "Switch" + str(label) + "%")
        
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

def plot_err_mean(result_dir, save_dir, switching_values, my_colors):
    figsize = (4, 3)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    
    MAEs = {}
    RMSEs = {}
    MAPEs = {}
    xlabels = []
    
    sum_err_S = {}
    sum_err_I = {}
    sum_err_R = {}
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/means.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
            sum_err_S[sv] = 0.
            sum_err_I[sv] = 0.
            sum_err_R[sv] = 0.
            
    
    counter = 0
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0]
        label = f"{float(sv) * 100:.4f}%"
        error_S = (df['muS_r0'] - true['muS_r0']).abs()
        error_I = (df['muI_r0'] - true['muI_r0']).abs()
        error_R = (df['muR_r0'] - true['muR_r0']).abs()
        ax_S.plot(time, np.log(error_S), label = label, color=my_colors[counter], linewidth=0.5)
        ax_I.plot(time, np.log(error_I), label = label, color=my_colors[counter], linewidth=0.5)
        ax_R.plot(time, np.log(error_R), label = label, color=my_colors[counter], linewidth=0.5)
        sum_err_S[sv] = error_S.sum()
        sum_err_I[sv] = error_I.sum()
        sum_err_R[sv] = error_R.sum()
        counter += 1
        xlabels.append(label)
        
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
    ax_S.set_ylabel("Log(err(S))")
    #ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Log(err(I))")
    #ax_I.set_yscale("log")
    #ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Log(err(R))")
    #ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_R.png", dpi=dpi)
    plt.close(fig_R)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    labels = list(MAEs.keys())
    x = range(len(labels))
    
    # Bar plot with summed errors
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_S.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_S.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_I.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_I.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_R.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_R.png", dpi=dpi)
    plt.close()
    
    # Grouped bar plot with total MAE and MSE
    
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
        
def plot_err_var(result_dir, save_dir, switching_values, my_colors):
    figsize = (4 ,3)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    
    MAEs = {}
    RMSEs = {}
    MAPEs = {}  # added mapes
    
    sum_err_S = {}
    sum_err_I = {}
    sum_err_R = {}
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
            
    counter = 0
    xlabels = []
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0]
        label = f"{float(sv) * 100:.4f}%"
        error_S = (df['M200'] - true['M200']).abs()
        error_I = (df['M020'] - true['M020']).abs()
        error_R = (df['M002'] - true['M002']).abs()
        ax_S.plot(time, np.log(error_S), label = label, color=my_colors[counter], linewidth=5)
        ax_I.plot(time, np.log(error_I), label = label, color=my_colors[counter], linewidth=0.5)
        ax_R.plot(time, np.log(error_R), label = label, color=my_colors[counter], linewidth=0.5)
        sum_err_S[sv] = error_S.sum()
        sum_err_I[sv] = error_I.sum()
        sum_err_R[sv] = error_R.sum()
        counter += 1
        xlabels.append(label)
        
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
    ax_S.set_ylabel("Log(err(VarS))")
    #ax_S.legend()
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_VarS.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Log(err(VarI))")
    #ax_I.set_yscale("log")
    #ax_I.legend()
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_VarI.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Log(err(VarR))")
    #ax_R.legend()
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_VarR.png", dpi=dpi)
    plt.close(fig_R)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    labels = list(MAEs.keys())
    x = range(len(labels))
    
    # Bar plot with summed errors
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_S.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_VarS.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_I.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_VarI.png", dpi=dpi)
    plt.close()
    
    fig, ax = plt.subplots(figsize=figsize)
    val = list(sum_err_R.values())
    ax.bar(x, val, color = my_colors[:len(x)])
    ax.set_xticks(x)
    ax.set_xticklabels(["" for _ in range(len(x))])
    ax.set_ylabel("Summed error")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(save_dir + "Sum_Err_VarR.png", dpi=dpi)
    plt.close()
    
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
    figsize = (4,2.5)
    mean_runtimes = []
    max_runtimes = []
    xlabels = []
    for sv in switching_values:
        runtime = pd.read_csv(result_dir + f"switch_value_{sv}/runtimes.csv")
        mean_runtimes.append(runtime['runtime'].mean())
        max_runtimes.append(runtime['runtime'].max())
        #xlabels.append(f"{float(sv) * 100:.4f}%")
        xlabels.append("")
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(range(len(xlabels)), mean_runtimes, color=colors['dark red'])
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale('log')
    ax.set_xticklabels(["" for _ in range(len(xlabels))])
    ax.set_title("Mean runtime")
    fig.tight_layout()
    fig.savefig(save_dir + "runtimes_mean.png", dpi=dpi)
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(range(len(xlabels)), max_runtimes, color=colors['dark red'])
    ax.set_ylabel("Runtime [s]")
    ax.set_yscale('log')
    ax.set_xticklabels(["" for _ in range(len(xlabels))])
    ax.set_title("Max runtime")
    fig.tight_layout()
    fig.savefig(save_dir + "runtimes_max.png", dpi=dpi)
    plt.close(fig)

def plot_err_mean_two_regions(result_dir, save_dir, switching_values, num_regions, my_colors):
    figsize = (3.3, 2.0)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    fig_S_rel, ax_S_rel = plt.subplots(figsize=figsize)
    fig_I_rel, ax_I_rel = plt.subplots(figsize=figsize)
    fig_R_rel, ax_R_rel = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    xlabels = []
    
    sum_err_S = {}
    sum_err_I = {}
    sum_err_R = {}
    
    sum_rel_err_S = {}
    sum_rel_err_I = {}
    sum_rel_err_R = {}
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/means.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
            sum_err_S[sv] = [0. for _ in range(num_regions)]
            sum_err_I[sv] = [0. for _ in range(num_regions)]
            sum_err_R[sv] = [0. for _ in range(num_regions)]
            
            sum_rel_err_S[sv] = [0. for _ in range(num_regions)]
            sum_rel_err_I[sv] = [0. for _ in range(num_regions)]
            sum_rel_err_R[sv] = [0. for _ in range(num_regions)]
            
    
    counter = 0
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0]
        label = f"{float(sv) * 100:.4f}%"
        for r in range(num_regions):
            error_S = (df[f'muS_r{r}'] - true[f'muS_r{r}']).abs()
            error_I = (df[f'muI_r{r}'] - true[f'muI_r{r}']).abs()
            error_R = (df[f'muR_r{r}'] - true[f'muR_r{r}']).abs()
            rel_error_S = ((df[f'muS_r{r}'] - true[f'muS_r{r}'])/true[f'muS_r{r}']).abs()
            rel_error_I = ((df[f'muI_r{r}'] - true[f'muI_r{r}'])/true[f'muI_r{r}']).abs()
            rel_error_R = ((df[f'muR_r{r}'] - true[f'muR_r{r}'])/true[f'muR_r{r}']).abs()
            ax_S.plot(time, np.log(error_S), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=5)
            ax_I.plot(time, np.log(error_I), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_R.plot(time, np.log(error_R), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_S_rel.plot(time, np.log(rel_error_S), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=5)
            ax_I_rel.plot(time, np.log(rel_error_I), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_R_rel.plot(time, np.log(rel_error_R), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            sum_err_S[sv][r] = error_S.sum()
            sum_err_I[sv][r] = error_I.sum()
            sum_err_R[sv][r] = error_R.sum()
            sum_rel_err_S[sv][r] = rel_error_S.sum()
            sum_rel_err_I[sv][r] = rel_error_I.sum()
            sum_rel_err_R[sv][r] = rel_error_R.sum()
            counter += 1
        xlabels.append(label)
        
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Log(err(S))")
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_S.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Log(err(I))")
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_I.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Log(err(R))")
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_R.png", dpi=dpi)
    plt.close(fig_R)
    
    ax_S_rel.set_xlabel("Time [days]")
    ax_S_rel.set_ylabel("Log(rel_err(S))")
    fig_S_rel.tight_layout()
    fig_S_rel.savefig(save_dir + "rel_Err_ts_S.png", dpi=dpi)
    plt.close(fig_S_rel)
    
    ax_I_rel.set_xlabel("Time [days]")
    ax_I_rel.set_ylabel("Log(rel_err(I))")
    fig_I_rel.tight_layout()
    fig_I_rel.savefig(save_dir + "rel_Err_ts_I.png", dpi=dpi)
    plt.close(fig_I_rel)
    
    ax_R_rel.set_xlabel("Time [days]")
    ax_R_rel.set_ylabel("Log(rel_err(R))")
    fig_R_rel.tight_layout()
    fig_R_rel.savefig(save_dir + "rel_Err_ts_R.png", dpi=dpi)
    plt.close(fig_R_rel)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    labels = list([f"{float(sv) * 100:.4f}%" for sv in dfs.keys()])
    x = range(len(labels))
    
    # Grouped bar plot    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=figsize)
    vals_r0_S = [sum_err_S[k][0] for k in sum_err_S.keys()]
    vals_r1_S = [sum_err_S[k][1] for k in sum_err_S.keys()]
    ax_bar_S.bar([i - 0.2 for i in x], vals_r0_S, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_S.bar([i + 0.2 for i in x], vals_r1_S, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_S.set_yscale('log')
    ax_bar_S.set_ylabel("Summed error")
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Sum_Err_S.png", dpi=dpi)
    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=figsize)
    vals_r0_S = [sum_rel_err_S[k][0] for k in sum_rel_err_S.keys()]
    vals_r1_S = [sum_rel_err_S[k][1] for k in sum_rel_err_S.keys()]
    ax_bar_S.bar([i - 0.2 for i in x], vals_r0_S, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_S.bar([i + 0.2 for i in x], vals_r1_S, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_S.set_yscale('log')
    ax_bar_S.set_ylabel("Summed rel error")
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Sum_rel_Err_S.png", dpi=dpi)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=figsize)
    vals_r0_I = [sum_err_I[k][0] for k in sum_err_I.keys()]
    vals_r1_I = [sum_err_I[k][1] for k in sum_err_I.keys()]
    ax_bar_I.bar([i - 0.2 for i in x], vals_r0_I, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_I.bar([i + 0.2 for i in x], vals_r1_I, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_I.set_yscale('log')
    ax_bar_I.set_ylabel("Summed error")
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Sum_Err_I.png", dpi=dpi)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=figsize)
    vals_r0_I = [sum_rel_err_I[k][0] for k in sum_rel_err_I.keys()]
    vals_r1_I = [sum_rel_err_I[k][1] for k in sum_rel_err_I.keys()]
    ax_bar_I.bar([i - 0.2 for i in x], vals_r0_I, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_I.bar([i + 0.2 for i in x], vals_r1_I, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_I.set_yscale('log')
    ax_bar_I.set_ylabel("Summed rel error")
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Sum_rel_Err_I.png", dpi=dpi)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=figsize)
    vals_r0_R = [sum_err_R[k][0] for k in sum_err_R.keys()]
    vals_r1_R = [sum_err_R[k][1] for k in sum_err_R.keys()]
    ax_bar_R.bar([i - 0.2 for i in x], vals_r0_R, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_R.bar([i + 0.2 for i in x], vals_r1_R, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_R.set_yscale('log')
    ax_bar_R.set_ylabel("Summed error")
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Sum_Err_R.png", dpi=dpi)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=figsize)
    vals_r0_R = [sum_rel_err_R[k][0] for k in sum_rel_err_R.keys()]
    vals_r1_R = [sum_rel_err_R[k][1] for k in sum_rel_err_R.keys()]
    ax_bar_R.bar([i - 0.2 for i in x], vals_r0_R, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_R.bar([i + 0.2 for i in x], vals_r1_R, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_R.set_yscale('log')
    ax_bar_R.set_ylabel("Summed rel error")
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Sum_rel_Err_R.png", dpi=dpi)

def plot_err_variance_two_regions(result_dir, save_dir, switching_values, num_regions, my_colors):
    figsize = (3.3, 2.0)
    
    fig_S, ax_S = plt.subplots(figsize=figsize)
    fig_I, ax_I = plt.subplots(figsize=figsize)
    fig_R, ax_R = plt.subplots(figsize=figsize)
    
    fig_S_rel, ax_S_rel = plt.subplots(figsize=figsize)
    fig_I_rel, ax_I_rel = plt.subplots(figsize=figsize)
    fig_R_rel, ax_R_rel = plt.subplots(figsize=figsize)
    
    dfs = {}
    true_df = {}
    xlabels = []
    
    sum_err_S = {}
    sum_err_I = {}
    sum_err_R = {}
    
    sum_rel_err_S = {}
    sum_rel_err_I = {}
    sum_rel_err_R = {}
    
    for sv in switching_values:
        means = pd.read_csv(result_dir + f"switch_value_{sv}/moments.csv")
        if(sv == "1.000000"):
            true_df[sv] = means
        else:
            dfs[sv] = means
            sum_err_S[sv] = [0. for _ in range(num_regions)]
            sum_err_I[sv] = [0. for _ in range(num_regions)]
            sum_err_R[sv] = [0. for _ in range(num_regions)]
            
            sum_rel_err_S[sv] = [0. for _ in range(num_regions)]
            sum_rel_err_I[sv] = [0. for _ in range(num_regions)]
            sum_rel_err_R[sv] = [0. for _ in range(num_regions)]
            
    
    counter = 0
    for sv, df in dfs.items():
        time = df["Time"]
        true = list(true_df.values())[0]
        label = f"{float(sv) * 100:.4f}%"
        for r in range(num_regions):
            string_S = "M"
            string_I = "M"
            string_R = "M"
            for _ in range(3*num_regions):
                string_S += "0"
                string_I += "0"
                string_R += "0"
            string_S = string_S[:(r * 3 + 1)] + "2" + string_S[(r * 3 + 1 + 1):]
            string_I = string_I[:(r * 3 + 2)] + "2" + string_I[(r * 3 + 2 + 1):]
            string_R = string_R[:(r * 3 + 3)] + "2" + string_R[(r * 3 + 3 + 1):]
            error_S = (df[string_S] - true[string_S]).abs()
            error_I = (df[string_I] - true[string_I]).abs()
            error_R = (df[string_R] - true[string_R]).abs()
            rel_error_S = ((df[string_S] - true[string_S])/true[string_S]).abs()
            rel_error_I = ((df[string_I] - true[string_I])/true[string_I]).abs()
            rel_error_R = ((df[string_R] - true[string_R])/true[string_R]).abs()
            ax_S.plot(time, np.log(error_S), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=5)
            ax_I.plot(time, np.log(error_I), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_R.plot(time, np.log(error_R), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_S_rel.plot(time, np.log(rel_error_S), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=5)
            ax_I_rel.plot(time, np.log(rel_error_I), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            ax_R_rel.plot(time, np.log(rel_error_R), label = str(label) + f"%, Region {r:.0f}", color=my_colors[counter], linewidth=0.5)
            sum_err_S[sv][r] = error_S.sum()
            sum_err_I[sv][r] = error_I.sum()
            sum_err_R[sv][r] = error_R.sum()
            sum_rel_err_S[sv][r] = rel_error_S.sum()
            sum_rel_err_I[sv][r] = rel_error_I.sum()
            sum_rel_err_R[sv][r] = rel_error_R.sum()
            counter += 1
        xlabels.append(label)
        
        
    ax_S.set_xlabel("Time [days]")
    ax_S.set_ylabel("Log(err(VarS))")
    fig_S.tight_layout()
    fig_S.savefig(save_dir + "Err_ts_VarS.png", dpi=dpi)
    plt.close(fig_S)
    
    ax_I.set_xlabel("Time [days]")
    ax_I.set_ylabel("Log(err(VarI))")
    fig_I.tight_layout()
    fig_I.savefig(save_dir + "Err_ts_VarI.png", dpi=dpi)
    plt.close(fig_I)
    
    ax_R.set_xlabel("Time [days]")
    ax_R.set_ylabel("Log(err(VarR))")
    fig_R.tight_layout()
    fig_R.savefig(save_dir + "Err_ts_VarR.png", dpi=dpi)
    plt.close(fig_R)
    
    ax_S_rel.set_xlabel("Time [days]")
    ax_S_rel.set_ylabel("Log(rel_err(VarS))")
    fig_S_rel.tight_layout()
    fig_S_rel.savefig(save_dir + "rel_Err_ts_VarS.png", dpi=dpi)
    plt.close(fig_S_rel)
    
    ax_I_rel.set_xlabel("Time [days]")
    ax_I_rel.set_ylabel("Log(rel_err(VarI))")
    fig_I_rel.tight_layout()
    fig_I_rel.savefig(save_dir + "rel_Err_ts_VarI.png", dpi=dpi)
    plt.close(fig_I_rel)
    
    ax_R_rel.set_xlabel("Time [days]")
    ax_R_rel.set_ylabel("Log(rel_err(VarR))")
    fig_R_rel.tight_layout()
    fig_R_rel.savefig(save_dir + "rel_Err_ts_VarR.png", dpi=dpi)
    plt.close(fig_R_rel)
    
    # Save legend seperately
    handles, labels = ax_S.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)                   
    fig_leg.legend(handles, labels, loc='center', ncol=2)      
    fig_leg.tight_layout()
    fig_leg.savefig(save_dir + "legend_S_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_I.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_I_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

    handles, labels = ax_R.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=figsize)
    fig_leg.legend(handles, labels, loc='center', ncol=2)
    fig_leg.savefig(save_dir + "legend_R_err.png", dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig_leg)
    
    labels = list([f"{float(sv) * 100:.4f}%" for sv in dfs.keys()])
    x = range(len(labels))
    
    # Grouped bar plot    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=figsize)
    vals_r0_S = [sum_err_S[k][0] for k in sum_err_S.keys()]
    vals_r1_S = [sum_err_S[k][1] for k in sum_err_S.keys()]
    ax_bar_S.bar([i - 0.2 for i in x], vals_r0_S, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_S.bar([i + 0.2 for i in x], vals_r1_S, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_S.set_yscale('log')
    ax_bar_S.set_ylabel("Summed error")
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Sum_Err_VarS.png", dpi=dpi)
    
    fig_bar_S, ax_bar_S = plt.subplots(figsize=figsize)
    vals_r0_S = [sum_rel_err_S[k][0] for k in sum_rel_err_S.keys()]
    vals_r1_S = [sum_rel_err_S[k][1] for k in sum_rel_err_S.keys()]
    ax_bar_S.bar([i - 0.2 for i in x], vals_r0_S, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_S.bar([i + 0.2 for i in x], vals_r1_S, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_S.set_xticks(x)
    ax_bar_S.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_S.set_yscale('log')
    ax_bar_S.set_ylabel("Summed rel error")
    fig_bar_S.tight_layout()
    fig_bar_S.savefig(save_dir + "Sum_rel_Err_VarS.png", dpi=dpi)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=figsize)
    vals_r0_I = [sum_err_I[k][0] for k in sum_err_I.keys()]
    vals_r1_I = [sum_err_I[k][1] for k in sum_err_I.keys()]
    ax_bar_I.bar([i - 0.2 for i in x], vals_r0_I, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_I.bar([i + 0.2 for i in x], vals_r1_I, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_I.set_yscale('log')
    ax_bar_I.set_ylabel("Summed error")
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Sum_Err_VarI.png", dpi=dpi)
    
    fig_bar_I, ax_bar_I = plt.subplots(figsize=figsize)
    vals_r0_I = [sum_rel_err_I[k][0] for k in sum_rel_err_I.keys()]
    vals_r1_I = [sum_rel_err_I[k][1] for k in sum_rel_err_I.keys()]
    ax_bar_I.bar([i - 0.2 for i in x], vals_r0_I, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_I.bar([i + 0.2 for i in x], vals_r1_I, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_I.set_xticks(x)
    ax_bar_I.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_I.set_yscale('log')
    ax_bar_I.set_ylabel("Summed rel error")
    fig_bar_I.tight_layout()
    fig_bar_I.savefig(save_dir + "Sum_rel_Err_VarI.png", dpi=dpi)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=figsize)
    vals_r0_R = [sum_err_R[k][0] for k in sum_err_R.keys()]
    vals_r1_R = [sum_err_R[k][1] for k in sum_err_R.keys()]
    ax_bar_R.bar([i - 0.2 for i in x], vals_r0_R, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_R.bar([i + 0.2 for i in x], vals_r1_R, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_R.set_yscale('log')
    ax_bar_R.set_ylabel("Summed error")
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Sum_Err_VarR.png", dpi=dpi)
    
    fig_bar_R, ax_bar_R = plt.subplots(figsize=figsize)
    vals_r0_R = [sum_rel_err_R[k][0] for k in sum_rel_err_R.keys()]
    vals_r1_R = [sum_rel_err_R[k][1] for k in sum_rel_err_R.keys()]
    ax_bar_R.bar([i - 0.2 for i in x], vals_r0_R, width=0.4, color = [my_colors[c] for c in range(0, len(my_colors), 2)])
    ax_bar_R.bar([i + 0.2 for i in x], vals_r1_R, width=0.4, color = [my_colors[c] for c in range(1, len(my_colors), 2)])
    ax_bar_R.set_xticks(x)
    ax_bar_R.set_xticklabels(["" for _ in range(len(labels))])
    ax_bar_R.set_yscale('log')
    ax_bar_R.set_ylabel("Summed rel error")
    fig_bar_R.tight_layout()
    fig_bar_R.savefig(save_dir + "Sum_rel_Err_VarR.png", dpi=dpi)

config = "config_2r1"         
result_dir = "V:/bick_ju/TemporalHybrid/Hybrid1/" + config + "/"
switching_values = ["0.000000", "0.001000" ,"0.010000", "0.100000", "1.000000"]
save_dir = "H:/Documents/TemporalHybridModel/Hybrid1/" + config + "/"
os.makedirs(save_dir, exist_ok=True)

plot_means(result_dir, save_dir, switching_values, 2, [colors['dark blue'], colors['purple'], colors['middle blue'], colors['rose'], colors['light blue'], colors['red'], colors['light teal'], colors['brown'], colors['dark green'], colors['middle green']])
plot_variances(result_dir, save_dir, switching_values, 2, [colors['dark blue'], colors['purple'], colors['middle blue'], colors['rose'], colors['light blue'], colors['red'], colors['light teal'], colors['brown'], colors['dark green'], colors['middle green']])

# plot_covariances(result_dir, save_dir, switching_values)
# plot_err_mean(result_dir, save_dir, switching_values, list(colors.values()))
# plot_err_var(result_dir, save_dir, switching_values, list(colors.values()))
# plot_runtimes(result_dir, save_dir, switching_values)

plot_err_mean_two_regions(result_dir, save_dir, switching_values, 2, [colors['dark blue'], colors['purple'], colors['middle blue'], colors['rose'], colors['light blue'], colors['red'], colors['light teal'], colors['brown'], colors['dark green'], colors['middle green']])
plot_err_variance_two_regions(result_dir, save_dir, switching_values, 2, [colors['dark blue'], colors['purple'], colors['middle blue'], colors['rose'], colors['light blue'], colors['red'], colors['light teal'], colors['brown'], colors['dark green'], colors['middle green']])
