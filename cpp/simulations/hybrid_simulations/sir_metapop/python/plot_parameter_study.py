from settings import *
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import LogLocator, NullLocator


def format_bin_edge(value):
    if value > 0 and float(value).is_integer():
        log = np.log10(value)
        if (np.isclose(log, round(log))) and (log != 0):
            return r"$10^{%d}$" % round(log)
        if (value % 10 == 0):
            return fr"${int(value/10**int(log))}\cdot 10^{int(log)}$"
    return f"{int(value)}" if float(value).is_integer() else f"{value}"


def bin_labels(edges):
    return [
        fr"{format_bin_edge(edges[1])}" if e == -np.inf else
        fr"$\geq${format_bin_edge(edges[-2])}" if edges[i + 1] == np.inf else
        f"{format_bin_edge(edges[i + 1])}"
        for i, e in enumerate(edges[:-1])
    ]


def plot_parameter_distributions(df, save_path, figsize):
    params = ["I_init_mean", "R_init_mean", "lambda",
              "R0_mean", "std_I_init", "std_R_init", "std_I_init_norm"]
    param_names = [r"$\mu_I(t_{init})$", r"$\mu_R(t_{init})$", r"$\lambda$",
                   r"$R_0(\mu(t_{init}))$", r"$\sigma_I(t_{init})$", r"$\sigma_R(t_{init})$", r"$\sigma_I(t_{init})/\mu_I(t_{init})$"]

    # per-parameter bin edges; use -np.inf / np.inf for open-ended bins
    bin_edges = {
        "I_init_mean": [-np.inf, 10, 100, 1000, 10000, np.inf],
        "R_init_mean": [-np.inf, 1000, 10000, 50000, np.inf],
        "lambda": 20,
        "R0_mean": [-np.inf, 1, 2, 3, np.inf],
        "std_I_init": [-np.inf, 10, 100, 1000, np.inf],
        "std_R_init": [-np.inf, 10, 100, 1000, np.inf],
        "std_I_init_norm": [-np.inf, 0, 0.7, 1, 2, np.inf]
    }

    for param, name in zip(params, param_names):
        fig_bins, ax_bins = plt.subplots(figsize=figsize)
        fig, ax = plt.subplots(figsize=figsize)
        edges = bin_edges[param]

        if isinstance(edges, int):
            ax_bins.hist(df[param].dropna(), bins=edges,
                         color=colors["purple"])
        else:
            labels = bin_labels(edges)
            counts = pd.cut(df[param].dropna(), bins=edges,
                            labels=labels).value_counts().reindex(labels)
            ax_bins.bar(labels, counts.values, color=colors["purple"])
            ax.hist(df[param].dropna(), bins=50, color=colors["purple"])

        ax_bins.set_title(name)
        ax_bins.set_ylabel("Frequency [#]")
        ax.set_title(name)
        ax.set_ylabel("Frequency [#]")
        formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((3, 3))
        ax_bins.yaxis.set_major_formatter(formatter)
        fig_bins.subplots_adjust(left=0.19, bottom=0.15, right=0.98)
        fig_bins.savefig(save_path + f"histogram_{param}.png", dpi=dpi)
        fig.tight_layout()
        fig.savefig(save_path + f"histogram_{param}_reg_bins.png", dpi=dpi)


def plot_parameter_nan_grid(df, save_path, figsize, quantity):
    params = ["I_init_mean", "R_init_mean", "R0_mean", "std_I_init_norm"]
    param_names = [r"$\mu_I(t_{init})$", r"$\mu_R(t_{init})$",
                   r"$R_0(\mu(t_{init}))$", r"$\sigma_I(t_{init})/\mu_I(t_{init})$"]
    n = len(params)

    # per-parameter bin edges; use -np.inf / np.inf for open-ended bins
    bin_edges = {
        "I_init_mean": [-np.inf, 10, 100, 1000, 10000, np.inf],
        "R_init_mean": [-np.inf, 1000, 10000, 50000, np.inf],
        "R0_mean": [-np.inf, 1, 2, 3, np.inf],
        "std_I_init_norm": [-np.inf, 0, 0.7, 1, 2, np.inf]
    }

    nan_weights = df[f"NaN_{quantity}"].astype(int)
    cmap = LinearSegmentedColormap.from_list(
        "nan_counts", ["white", colors["teal"], colors["dark teal"]])

    def custom_bins_1d(param):
        labels = bin_labels(bin_edges[param])
        cats = pd.cut(df[param], bins=bin_edges[param], labels=labels)
        counts = nan_weights.groupby(
            cats, observed=False).sum().reindex(labels)
        return labels, counts.values

    def custom_bins_2d(param_x, param_y):
        labels_x = bin_labels(bin_edges[param_x])
        labels_y = bin_labels(bin_edges[param_y])
        cats_x = pd.cut(df[param_x], bins=bin_edges[param_x], labels=labels_x)
        cats_y = pd.cut(df[param_y], bins=bin_edges[param_y], labels=labels_y)
        grid = (nan_weights.groupby([cats_y, cats_x], observed=False)
                .sum().unstack().reindex(index=labels_y, columns=labels_x))
        return labels_x, labels_y, grid.values

    def valid_x(param_x):
        mask = df[param_x].notna()
        return df[param_x][mask], nan_weights[mask]

    def valid_xy(param_x, param_y):
        mask = df[param_x].notna() & df[param_y].notna()
        return df[param_x][mask], df[param_y][mask], nan_weights[mask]

    for use_custom_edges, suffix in [(False, "reg_bins"), (True, "bin_edges")]:
        # shared color scale across all off-diagonal panels of this figure
        max_count = 0
        for i, param_x in enumerate(params):
            for j, param_y in enumerate(params):
                if i == j:
                    continue
                if use_custom_edges:
                    *_, grid = custom_bins_2d(param_x, param_y)
                else:
                    x, y, w = valid_xy(param_x, param_y)
                    grid, _, _ = np.histogram2d(x, y, bins=30, weights=w)
                    grid = np.log1p(grid)
                max_count = max(max_count, np.nanmax(grid))

        # shared y-scale across all diagonal panels of this figure
        max_diag = 0
        for param in params:
            if use_custom_edges:
                _, counts = custom_bins_1d(param)
            else:
                x, w = valid_x(param)
                counts, _ = np.histogram(x, bins=20, weights=w)
            max_diag = max(max_diag, np.nanmax(counts))

        fig, axes = plt.subplots(
            n, n, figsize=(figsize[0] * n, figsize[1] * n))
        im = None
        for i, param_x in enumerate(params):
            for j, param_y in enumerate(params):
                ax = axes[j, i]

                if i == j:
                    if use_custom_edges:
                        labels, counts = custom_bins_1d(param_x)
                        ax.bar(labels, counts, color=colors["middle blue"])
                        ax.tick_params(axis="x", rotation=45)
                        ax.set_ylim(0, max_diag+10)
                        if (i != 0):
                            ax.set_yticklabels([])
                        if (i != n-1):
                            ax.set_xticklabels([])
                    else:
                        x, w = valid_x(param_x)
                        ax.hist(x, bins=20, weights=w,
                                color=colors["middle blue"])
                        ax.set_ylim(0, max_diag)
                else:
                    if use_custom_edges:
                        labels_x, labels_y, grid = custom_bins_2d(
                            param_x, param_y)
                        im = ax.pcolormesh(
                            grid, cmap=cmap, vmin=0, vmax=max_count)
                        ax.set_xticks(np.arange(len(labels_x)) + 0.5)
                        ax.set_yticks(np.arange(len(labels_y)) + 0.5)
                        ax.set_xticklabels(
                            labels_x if j == n - 1 else [],
                            rotation=45, ha="right")
                        ax.set_yticklabels(labels_y if i == 0 else [])
                    else:
                        x, y, w = valid_xy(param_x, param_y)
                        grid, xedges, yedges = np.histogram2d(
                            x, y, bins=20, weights=w)
                        im = ax.pcolormesh(
                            xedges, yedges, np.log1p(grid).T,
                            cmap=cmap, vmin=0, vmax=max_count)
                        if j != n - 1:
                            ax.set_xticklabels([])
                        if i != 0:
                            ax.set_yticklabels([])

                if j == n - 1:
                    ax.set_xlabel(param_names[i])
                if i == 0:
                    ax.set_ylabel(param_names[j])

        assert im is not None
        cbar_label = "NaN values [#]" if use_custom_edges else "log(NaN count + 1)"
        cbar = fig.colorbar(im, ax=axes, fraction=0.03,
                            pad=0.5, shrink=0.8, anchor=(2, 0.7))
        cbar.set_label(cbar_label, labelpad=15)
        fig.subplots_adjust(bottom=0.13, right=0.85, left=0.09, top=0.97)
        fig.savefig(save_path + f"nan_grid_{suffix}_{quantity}.png", dpi=dpi)


def plot_parameter_err_grid(df, save_path, figsize, quantity, error, log_scale=False):
    params = ["I_init_mean", "R_init_mean", "R0_mean", "std_I_init_norm"]
    param_names = [r"$\mu_I(t_{init})$", r"$\mu_R(t_{init})$",
                   r"$R_0(\mu(t_{init}))$", r"$\sigma_I(t_{init})/\mu_I(t_{init})$"]
    n = len(params)

    # per-parameter bin edges; use -np.inf / np.inf for open-ended bins
    bin_edges = {
        "I_init_mean": [-np.inf, 10, 100, 1000, 10000, np.inf],
        "R_init_mean": [-np.inf, 1000, 10000, 50000, np.inf],
        "R0_mean": [-np.inf, 1, 2, 3, np.inf],
        "std_I_init_norm": [-np.inf, 0, 0.7, 1, 2, np.inf]
    }

    smape = df[f"{error}_{quantity}"]
    cmap = LinearSegmentedColormap.from_list(
        "smape", [colors["light teal"], colors["teal"], colors["dark teal"]])

    def custom_bins_1d(param):
        labels = bin_labels(bin_edges[param])
        cats = pd.cut(df[param], bins=bin_edges[param], labels=labels)
        means = smape.groupby(cats, observed=False).mean().reindex(labels)
        return labels, means.values

    def custom_bins_2d(param_x, param_y):
        labels_x = bin_labels(bin_edges[param_x])
        labels_y = bin_labels(bin_edges[param_y])
        cats_x = pd.cut(df[param_x], bins=bin_edges[param_x], labels=labels_x)
        cats_y = pd.cut(df[param_y], bins=bin_edges[param_y], labels=labels_y)
        grid = (smape.groupby([cats_y, cats_x], observed=False)
                .mean().unstack().reindex(index=labels_y, columns=labels_x))
        return labels_x, labels_y, grid.values

    def valid_x(param_x):
        mask = df[param_x].notna() & smape.notna()
        return df[param_x][mask], smape[mask]

    def valid_xy(param_x, param_y):
        mask = df[param_x].notna() & df[param_y].notna() & smape.notna()
        return df[param_x][mask], df[param_y][mask], smape[mask]

    def binned_mean_1d(x, w, bins=20):
        counts, edges = np.histogram(x, bins=bins)
        sums, _ = np.histogram(x, bins=edges, weights=w)
        with np.errstate(invalid="ignore"):
            means = np.where(counts > 0, sums / counts, np.nan)
        return means, edges

    def binned_mean_2d(x, y, w, bins=20):
        counts, xedges, yedges = np.histogram2d(x, y, bins=bins)
        sums, _, _ = np.histogram2d(
            x, y, bins=[xedges, yedges], weights=w)
        with np.errstate(invalid="ignore"):
            means = np.where(counts > 0, sums / counts, np.nan)
        return means, xedges, yedges

    for use_custom_edges, suffix in [(False, "reg_bins"), (True, "bin_edges")]:
        # shared color scale across all off-diagonal panels of this figure
        max_smape = 0
        min_positive_smape = np.inf
        for i, param_x in enumerate(params):
            for j, param_y in enumerate(params):
                if i == j:
                    continue
                if use_custom_edges:
                    *_, grid = custom_bins_2d(param_x, param_y)
                else:
                    x, y, w = valid_xy(param_x, param_y)
                    grid, _, _ = binned_mean_2d(x, y, w, bins=20)
                max_smape = max(max_smape, np.nanmax(grid))
                positive = grid[grid > 0]
                if positive.size:
                    min_positive_smape = min(
                        min_positive_smape, np.nanmin(positive))

        if log_scale:
            norm = LogNorm(vmin=min_positive_smape, vmax=max_smape)
        else:
            norm = Normalize(vmin=0, vmax=max_smape)

        # shared y-scale across all diagonal panels of this figure
        max_diag = 0
        min_positive_diag = np.inf
        for param in params:
            if use_custom_edges:
                _, means = custom_bins_1d(param)
            else:
                x, w = valid_x(param)
                means, _ = binned_mean_1d(x, w, bins=20)
            max_diag = max(max_diag, np.nanmax(means))
            positive = means[means > 0]
            if positive.size:
                min_positive_diag = min(
                    min_positive_diag, np.nanmin(positive))

        fig, axes = plt.subplots(
            n, n, figsize=(figsize[0] * n, figsize[1] * n))
        im = None
        for i, param_x in enumerate(params):
            for j, param_y in enumerate(params):
                ax = axes[j, i]

                if i == j:
                    ylim = (min_positive_diag, max_diag*1.2) if log_scale \
                        else (0, max_diag*1.01)
                    if log_scale:
                        ax.set_yscale("log")
                        ax.yaxis.set_major_locator(LogLocator(base=10.0))
                        ax.yaxis.set_minor_locator(NullLocator())
                    if use_custom_edges:
                        labels, means = custom_bins_1d(param_x)
                        ax.bar(labels, means, color=colors["middle blue"])
                        ax.tick_params(axis="x", rotation=45)
                        ax.set_ylim(*ylim)
                        if (i != 0):
                            ax.set_yticklabels([])
                        if (i != n-1):
                            ax.set_xticklabels([])
                    else:
                        x, w = valid_x(param_x)
                        means, edges = binned_mean_1d(x, w, bins=20)
                        ax.stairs(means, edges, fill=True,
                                  color=colors["middle blue"])
                        ax.set_ylim(*ylim)
                else:
                    if use_custom_edges:
                        labels_x, labels_y, grid = custom_bins_2d(
                            param_x, param_y)
                        im = ax.pcolormesh(
                            grid, cmap=cmap, norm=norm)
                        ax.set_xticks(np.arange(len(labels_x)) + 0.5)
                        ax.set_yticks(np.arange(len(labels_y)) + 0.5)
                        ax.set_xticklabels(
                            labels_x if j == n - 1 else [],
                            rotation=45, ha="right")
                        ax.set_yticklabels(labels_y if i == 0 else [])
                    else:
                        x, y, w = valid_xy(param_x, param_y)
                        grid, xedges, yedges = binned_mean_2d(
                            x, y, w, bins=20)
                        im = ax.pcolormesh(
                            xedges, yedges, grid.T,
                            cmap=cmap, norm=norm)
                        if j != n - 1:
                            ax.set_xticklabels([])
                        if i != 0:
                            ax.set_yticklabels([])

                if j == n - 1:
                    ax.set_xlabel(param_names[i])
                if i == 0:
                    ax.set_ylabel(param_names[j])

        assert im is not None
        cbar = fig.colorbar(im, ax=axes, fraction=0.03,
                            pad=0.5, shrink=0.8, anchor=(2, 0.7))
        err_quantity = ""
        if quantity == "mean":
            err_quantity = r"$(\mu_I)$"
        elif quantity == "var":
            err_quantity = r"$(\sigma^2_I)$"

        cbar.set_label(f"Mean {error}{err_quantity}", labelpad=10)
        log_suffix = "_log" if log_scale else ""
        fig.subplots_adjust(bottom=0.13, right=0.85, left=0.09, top=0.97)
        fig.savefig(
            save_path + f"{error}_grid_{suffix}_{quantity}{log_suffix}.png", dpi=dpi)


base_folder = "/Users/julia/sim_outputs/output/ParameterStudy/SIR/"
df = pd.read_csv(base_folder + "results.csv")
figsize = (4, 2.5)

plot_parameter_distributions(df, base_folder, figsize=(2.7, 2))
# plot_parameter_nan_grid(df, base_folder, figsize=figsize, quantity="mean")
# plot_parameter_nan_grid(df, base_folder, figsize=figsize, quantity="var")
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="mean", error="sMAPE")
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="var", error="sMAPE")
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="mean", error="MSE")
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="var", error="MSE")
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="mean", error="MSE", log_scale=True)
# plot_parameter_err_grid(df, base_folder, figsize=figsize,
#                         quantity="var", error="MSE", log_scale=True)
