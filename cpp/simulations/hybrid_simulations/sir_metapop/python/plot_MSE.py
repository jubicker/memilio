import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
import pandas as pd
import seaborn as sns


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


def plot_err_2D_v2(df, col, param1, param2):
    # Pivot to create heatmap grid
    heatmap_data = df.pivot_table(
        index=param1,
        columns=param2,
        values=col
    )

    # Mask NaNs for heatmap coloring
    mask = heatmap_data.isna()

    plt.figure(figsize=(10, 8))

    # Heatmap
    ax = sns.heatmap(
        heatmap_data,
        mask=mask,
        cmap='viridis',
        cbar_kws={'label': col}
    )

    # --- Overlay NaN points as red dots ---
    nan_df = df[df[col].isna() | df[param1].isna() | df[param2].isna()]

    # IMPORTANT: heatmap uses categorical grid positions, so we map values to indices
    a_order = heatmap_data.index.tolist()
    b_order = heatmap_data.columns.tolist()

    a_to_i = {v: i + 2 for i, v in enumerate(a_order)}
    b_to_j = {v: j + 2 for j, v in enumerate(b_order)}

    x = nan_df[param2].map(b_to_j)
    y = nan_df[param1].map(a_to_i)

    ax.scatter(x, y, color='red', s=50, label='NaN (invalid)')

    plt.xlabel(param1)
    plt.ylabel(param2)
    plt.legend()

    plt.savefig(dir + col + "_" + param1 +
                "_" + param2 + "_v2" + ".png")


dir = "/Users/julia/sim_outputs/output/MSE/SIR/"

mean_or_var = "var"
col = "sMAPE_" + mean_or_var

df = pd.read_csv(dir + "MSEs.csv")

plot_MSE_3D(df, col)
plot_MSE_3D(df, col, "max_relation", "max_mean", "R_0")

plot_MSE_2D(df, col, "max_relation_adapted", "R_0")

plot_err_2D_v2(df, col, "max_relation", "R_0")
