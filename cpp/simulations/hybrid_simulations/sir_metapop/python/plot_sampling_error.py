from settings import *
import os


def plot_error_histograms(err_file, err_file1, save_dir, title, labels=("Naive", "Rounding correction"), bins=10, figsize=(10, 4)):
    err = pd.read_csv(err_file).iloc[:, 0]
    err1 = pd.read_csv(err_file1).iloc[:, 0]

    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    axes[0].hist(err, bins=bins, color=colors["middle blue"], label=labels[0],
                 weights=100*np.ones_like(err) / len(err))
    mean = err.mean()
    axes[0].axvline(mean, color=colors["dark blue"], linestyle="--",
                    label=rf"$\mu_{{\text{{err}}}}^{{(\text{{naive}})}}$ = {mean:.4f}")
    axes[0].set_ylabel("Relative\nfrequency [%]", labelpad=13)
    axes[0].legend()

    axes[1].hist(err1, bins=bins, color=colors["middle green"], label=labels[1],
                 weights=100*np.ones_like(err1) / len(err1))
    mean1 = err1.mean()
    axes[1].axvline(mean1, color=colors["dark green"], linestyle="--",
                    label=rf"$\mu_{{\text{{err}}}}^{{(\text{{corrected}})}}$ = {mean1:.4f}")
    axes[1].set_xlabel("Error")
    axes[1].set_ylabel("Relative\nfrequency [%]", labelpad=13)
    axes[1].legend()

    fig.suptitle(title)
    fig.subplots_adjust(hspace=0.3, bottom=0.1,
                        top=0.92, left=0.24, right=0.95)
    fig.savefig(f"{save_dir}/err_histograms_{title}.png", dpi=dpi)
    plt.close(fig)


dir = "/Users/julia/repos/fork/memilio/output"
save_dir = "/Users/julia/repos/fork/memilio/output"
hybrid_model = "Spatial-Hybrid2"
config = "config_SIR_R0_1_1.5_2_4_exchange"
condition = "combined_relation_var_gradient_condition_region"
closure_method = "truncation"
closure_order = "closure_order_3"
smm_dir = f"{dir}/SMM/{config}"
hybrid_dir = f"{dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
save_dir = f"{save_dir}/{hybrid_model}/{config}/{condition}/{closure_method}/{closure_order}"
os.makedirs(save_dir, exist_ok=True)

err_file_var = f"{hybrid_dir}/err_var.csv"
err_file_corr_var = f"{hybrid_dir}/err_var_rounding_correction.csv"

err_file_mean = f"{hybrid_dir}/err_mean.csv"
err_file_corr_mean = f"{hybrid_dir}/err_mean_rounding_correction.csv"


plot_error_histograms(err_file_var, err_file_corr_var, save_dir,
                      title="Variance", figsize=(3.8, 5))
plot_error_histograms(err_file_mean, err_file_corr_mean, save_dir,
                      title="Mean", figsize=(3.8, 5))
