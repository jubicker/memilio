import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter


def plot_MSE_3D(df, col, param1="I_init", param2="R_init", param3="lambda"):

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Separate NaN and non-NaN values
    valid = df[col].notna()
    nan_mask = df[col].isna()

    vmin = df.loc[valid, col].min()
    vmax = df.loc[valid, col].max()

    if vmin <= 0:
        norm = SymLogNorm(
            linthresh=1,   # linear range around 0
            linscale=1,
            vmin=vmin,
            vmax=vmax
        )
    else:
        norm = LogNorm(
            vmin=df.loc[valid, col].min(),
            vmax=df.loc[valid, col].max()
        )

    # Plot valid points with color according to MSE
    sc = ax.scatter(
        df.loc[valid, param1],
        df.loc[valid, param2],
        df.loc[valid, param3],
        c=df.loc[valid, col],
        cmap='viridis',
        norm=norm
    )

    # Plot NaN points in white
    ax.scatter(
        df.loc[nan_mask, param1],
        df.loc[nan_mask, param2],
        df.loc[nan_mask, param3],
        color='red',
        edgecolor='red'
    )

    # Move colorbar further right
    cbar = plt.colorbar(sc, ax=ax, pad=0.15)
    cbar.set_label(col)

    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_zlabel(param3)

    plt.savefig(dir + col + "_" + param1 +
                "_" + param2 + "_" + param3 + ".png")


def plot_MSE_2D(df, col, param1="I_init", param2="R_init"):

    fig, ax = plt.subplots(figsize=(8, 6))

    # Separate NaN and non-NaN values
    valid = df[col].notna()
    nan_mask = df[col].isna()

    vmin = df.loc[valid, col].min()
    vmax = df.loc[valid, col].max()

    if vmin <= 0:
        norm = SymLogNorm(
            linthresh=1,   # linear range around 0
            linscale=1,
            vmin=vmin,
            vmax=vmax
        )
    else:
        norm = LogNorm(
            vmin=df.loc[valid, col].min(),
            vmax=df.loc[valid, col].max()
        )

    # Plot valid points with color according to MSE
    sc = ax.scatter(
        df.loc[valid, param1],
        df.loc[valid, param2],
        c=df.loc[valid, col],
        cmap='viridis',
        norm=norm
    )

    # Plot NaN points in white
    ax.scatter(
        df.loc[nan_mask, param1],
        df.loc[nan_mask, param2],
        color='red',
        edgecolor='red'
    )

    # Move colorbar further right
    cbar = plt.colorbar(sc, ax=ax, pad=0.15)
    cbar.set_label(col)

    ax.set_xlabel(param1)
    ax.set_ylabel(param2)

    plt.savefig(dir + col + "_" + param1 +
                "_" + param2 + "_" + ".png")


def plot_err_2D_bins(df, col, param1, param2, bins1, bins2, bin1_labels, bin2_labels):
    df[param1 + "_bin"] = pd.cut(
        df[param1],
        bins=bins1,
        labels=bin1_labels,
        right=False
    )

    df[param2 + "_bin"] = pd.cut(
        df[param2],
        bins=bins2,
        labels=bin2_labels,
        right=False
    )
    # -------------------------
    # Heatmap 1: Mean Err
    # -------------------------

    valid_df = df.dropna(subset=[col])

    mean_err = (
        valid_df
        .groupby([param2 + "_bin", param1 + "_bin"], observed=False)[col]
        .mean()
        .unstack()
    )

    # -------------------------
    # Heatmap 2: Number of NaNs
    # -------------------------

    nan_counts = (
        df.assign(is_nan=df[col].isna())
        .groupby([param2 + "_bin", param1 + "_bin"], observed=False)["is_nan"]
        .sum()
        .unstack()
    )

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(
        mean_err,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        ax=ax
    )

    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_title("Average " + col + " by parameter bins")
    plt.savefig(dir + col + "_" + param1 + "_" + param2 + "_bins_mean_err.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(
        nan_counts,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        ax=ax
    )

    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_title("Number of NaNs in " + col + " by parameter bins")
    plt.savefig(dir + col + "_" + param1 + "_" +
                param2 + "_bins_nan_counts.png")
    plt.close(fig)


def plot_err_3D_bins(df, col, param1, param2, param3, bins1, bins2, bins3, bin1_labels, bin2_labels, bin3_labels, z_coordinate_scientific):

    df[param1 + "_bin"] = pd.cut(
        df[param1],
        bins=bins1,
        labels=bin1_labels,
        right=False
    )

    df[param2 + "_bin"] = pd.cut(
        df[param2],
        bins=bins2,
        labels=bin2_labels,
        right=False
    )

    df[param3 + "_bin"] = pd.cut(
        df[param3],
        bins=bins3,
        labels=bin3_labels,
        right=False
    )

    grouped = (
        df.dropna(subset=[col])
        .groupby(
            [param1 + "_bin", param2 + "_bin", param3 + "_bin"],
            observed=False
        )[col]
        .mean()
        .reset_index()
    )

    x = grouped[param1 + "_bin"].cat.codes
    y = grouped[param2 + "_bin"].cat.codes
    z = grouped[param3 + "_bin"].cat.codes

    err = grouped[col]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    dx = dy = dz = 0.5

    norm = plt.Normalize(err.min(), err.max())
    colors = plt.cm.viridis(norm(err))

    ax.bar3d(x, y, z, dx, dy, dz, color=colors)

    mappable = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    plt.colorbar(mappable, ax=ax, label="Average " + col)

    ax.set_xlabel(param1, labelpad=10)
    ax.set_ylabel(param2)
    ax.set_zlabel(param3)

    ax.set_xticks(range(len(bin1_labels)))
    ax.set_xticklabels(bin1_labels, rotation=45, ha="right")

    ax.set_yticklabels([])
    for i, label in enumerate(bin2_labels):
        ax.text(
            6.25, i - 0.3, 0,
            str(label),
            zdir=None,
        )
    ax.set_zticks(range(len(bin3_labels)))
    ax.set_zticklabels(bin3_labels)
    if z_coordinate_scientific:
        ax.text(6.25, 6.25, 7, r"$\times 10^{-6}$", zdir=None)

    plt.savefig(dir + col + "_" + param1 + "_" +
                param2 + "_" + param3 + "_bins_3D.png")
    plt.close(fig)

    # Nan counts for 3D bins
    grouped = (
        df.assign(is_nan=df[col].isna())
        .groupby(
            [param1 + "_bin", param2 + "_bin", param3 + "_bin"],
            observed=False
        )["is_nan"]
        .sum()
        .reset_index()
    )
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    for c in [param1 + "_bin", param2 + "_bin", param3 + "_bin"]:
        grouped[c] = pd.Categorical(grouped[c], ordered=True)

    # Select only bins with non-zero counts for plotting
    plot_df = grouped

    x_plot = plot_df[param1 + "_bin"].cat.codes.values
    y_plot = plot_df[param2 + "_bin"].cat.codes.values
    z_plot = plot_df[param3 + "_bin"].cat.codes.values

    counts_plot = plot_df["is_nan"].values

    norm = plt.Normalize(counts_plot.min(), counts_plot.max())
    colors = plt.cm.viridis(norm(counts_plot))
    # make zeros fully transparent
    alpha = np.where(counts_plot == 0, 0.0, 1.0)

    colors = np.concatenate([colors[:, :3], alpha[:, None]], axis=1)

    ax.bar3d(x_plot, y_plot, z_plot, dx, dy, dz, color=colors)

    mappable = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    plt.colorbar(mappable, ax=ax, label="NaN counts " + col)

    ax.set_xlabel(param1, labelpad=10)
    ax.set_ylabel(param2)
    ax.set_zlabel(param3)

    ax.set_xticks(range(len(bin1_labels)))
    ax.set_xticklabels(bin1_labels, rotation=45, ha="right")

    ax.set_yticklabels([])
    for i, label in enumerate(bin2_labels):
        ax.text(
            6.25, i - 0.3, 0,
            str(label),
            zdir=None,
        )
    ax.set_zticks(range(len(bin3_labels)))
    ax.set_zticklabels(bin3_labels)
    if z_coordinate_scientific:
        ax.text(6.25, 6.25, 7, r"$\times 10^{-6}$", zdir=None)

    fig.savefig(dir + col + "_" + param1 + "_" +
                param2 + "_" + param3 + "_bins_nan_3D.png")
    plt.close(fig)


def plot_err_bar_bins(df, col, param1, bins1, bin1_labels):

    df[param1 + "_bin"] = pd.cut(
        df[param1],
        bins=bins1,
        labels=bin1_labels,
        right=False
    )

    grouped = (
        df.dropna(subset=[col])
        .groupby(
            [param1 + "_bin"],
            observed=False
        )[col]
        .mean()
        .reset_index()
    )

    x = grouped[param1 + "_bin"].cat.codes

    err = grouped[col]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot()

    ax.bar(x, err)

    ax.set_xlabel(param1, labelpad=10)

    ax.set_xticks(range(len(bin1_labels)))
    ax.set_xticklabels(bin1_labels, rotation=45, ha="right")

    ax.set_ylabel("Average " + col)

    plt.savefig(dir + col + "_" + param1 + "_bins_bar.png")
    plt.close(fig)

    # Nan counts for 3D bins
    grouped = (
        df.assign(is_nan=df[col].isna())
        .groupby(
            [param1 + "_bin"],
            observed=False
        )["is_nan"]
        .sum()
        .reset_index()
    )
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot()

    for c in [param1 + "_bin"]:
        grouped[c] = pd.Categorical(grouped[c], ordered=True)

    # Select only bins with non-zero counts for plotting
    plot_df = grouped

    x_plot = plot_df[param1 + "_bin"].cat.codes.values

    counts_plot = plot_df["is_nan"].values

    ax.bar(x_plot, counts_plot)

    ax.set_xlabel(param1, labelpad=10)

    ax.set_xticks(range(len(bin1_labels)))
    ax.set_xticklabels(bin1_labels, rotation=45, ha="right")

    ax.set_ylabel("NaN [#]")

    fig.savefig(dir + col + "_" + param1 + "_bins_nan_bar.png")
    plt.close(fig)


dir = "V:/bick_ju/TemporalHybrid/MSE/SIR/"

mean_or_var = "var"
col = "sMAPE_" + mean_or_var

df = pd.read_csv(dir + "MSEs.csv")

# plot_MSE_3D(df, col)
# plot_MSE_3D(df, col, "max_relation", "relation_at_max_mu", "R_0")
# plot_MSE_3D(df, col, "relation_at_max_mu", "max_mean", "R_0")
# plot_MSE_3D(df, col, "max_relation", "max_mean", "R_0")

# plot_MSE_2D(df, col, "relation_at_max_mu", "R_0")
# plot_MSE_2D(df, col, "relation_at_max_mu", "max_relation")

parameters = ["I_init", "R_init", "lambda", "R_0"]
parameter_bins = {
    "I_init": [0, 10, 20, 30, 40, 50, np.inf],
    "R_init": [0, 0.01*100000, 0.1*100000, 0.2*100000, 0.3*100000, 0.4*100000, 0.5*100000, np.inf],
    "lambda": np.linspace(0.0000014, 0.000006, 7).tolist() + [np.inf],
    "R_0": [0, 1, 2, 3, 4, np.inf],
    "max_mean": [0, 0.01*100000, 0.1*100000, 0.2*100000, 0.3*100000, np.inf],
    "max_std": [0, 10, 100, 500, 1000, 3000],
    "max_std_normalized": [0, 0.1, 0.2, 0.3, 0.4, 0.5, np.inf],
}

bin_labels = {
    "I_init": [parameter_bins["I_init"][i] for i in range(len(parameter_bins["I_init"]) - 1)],
    "R_init": [int(parameter_bins["R_init"][i]/1) for i in range(len(parameter_bins["R_init"]) - 1)],
    "lambda": [np.round(parameter_bins["lambda"][i]*1e6, 1) for i in range(len(parameter_bins["lambda"]) - 1)],
    "R_0": [parameter_bins["R_0"][i] for i in range(len(parameter_bins["R_0"]) - 1)],
    "max_mean": [int(parameter_bins["max_mean"][i]/1) for i in range(len(parameter_bins["max_mean"]) - 1)],
    "max_std": [parameter_bins["max_std"][i] for i in range(len(parameter_bins["max_std"]) - 1)],
    "max_std_normalized": [parameter_bins["max_std_normalized"][i] for i in range(len(parameter_bins["max_std_normalized"]) - 1)]
}


plot_err_3D_bins(df, col, "I_init", "R_init", "lambda",
                 parameter_bins["I_init"], parameter_bins["R_init"], parameter_bins["lambda"],
                 bin_labels["I_init"], bin_labels["R_init"], bin_labels["lambda"], True)
plot_err_3D_bins(df, col, "I_init", "R_init", "R_0",
                 parameter_bins["I_init"], parameter_bins["R_init"], parameter_bins["R_0"],
                 bin_labels["I_init"], bin_labels["R_init"], bin_labels["R_0"], False)
plot_err_3D_bins(df, col, "I_init", "max_std", "max_std_normalized",
                 parameter_bins["I_init"], parameter_bins["max_std"], parameter_bins["max_std_normalized"],
                 bin_labels["I_init"], bin_labels["max_std"], bin_labels["max_std_normalized"], False)
plot_err_3D_bins(df, col, "max_mean", "max_std", "max_std_normalized",
                 parameter_bins["max_mean"], parameter_bins["max_std"], parameter_bins["max_std_normalized"],
                 bin_labels["max_mean"], bin_labels["max_std"], bin_labels["max_std_normalized"], False)
plot_err_bar_bins(df, col, "I_init",
                  parameter_bins["I_init"], bin_labels["I_init"])
plot_err_bar_bins(df, col, "R_init",
                  parameter_bins["R_init"], bin_labels["R_init"])
plot_err_bar_bins(df, col, "lambda",
                  parameter_bins["lambda"], bin_labels["lambda"])
plot_err_bar_bins(df, col, "R_0", parameter_bins["R_0"], bin_labels["R_0"])
plot_err_bar_bins(df, col, "max_mean",
                  parameter_bins["max_mean"], bin_labels["max_mean"])
plot_err_bar_bins(df, col, "max_std",
                  parameter_bins["max_std"], bin_labels["max_std"])
plot_err_bar_bins(df, col, "max_std_normalized",
                  parameter_bins["max_std_normalized"], bin_labels["max_std_normalized"])

# plot_MSE_3D(df, col, "max_mean", "max_std", "R_0")
# plot_MSE_3D(df, col, "max_mean", "max_std_normalized", "R_0")

# for param1 in parameters:
#     for param2 in parameters:
#         if param1 != param2:
#             bins1 = parameter_bins[param1]
#             bins2 = parameter_bins[param2]
#             plot_err_2D_bins(df, col, param1, param2, bins1, bins2, bins1[1:], bins2[1:])
