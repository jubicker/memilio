from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
from settings import *

save_dir = "H:/Documents/TemporalHybridModel"
os.makedirs(save_dir, exist_ok=True)

num_cores = [1, 2, 4, 8, 16, 32, 56]

# These are the runtimes for 5mio agents per region
# time_spatial_hybrid = [8453.97, 4400.92, 2361.28, 1246.44, 702.725, 409.678, 279.84]
# time_smm = [49868.9, 24767, 12919.8, 6478.65, 3551.08, 1809.49, 1064.94]

# These are the runtimes for 10mio agents per region
time_spatial_hybrid = [24499.4, 12773.6, 6765.56, 3448.72, 1867.35, 981.646, 590.3]
time_smm = [105455, 54994.3, 28639.9, 14397.7, 7662.86, 4004.42, 2280.91]

figsize = (6, 4)

speed_up_spatial_hybrid = [time_smm[i]/time_spatial_hybrid[i] for i in range(len(num_cores))]

fig, ax = plt.subplots(figsize=figsize)
fig2, ax2 = plt.subplots(figsize=figsize)
ax2.plot(num_cores, [time_smm[i]/60. for i in range(len(num_cores))], marker='o', color=colors["dark grey"], label="Stochastic metapopulation model")
ax2.plot(num_cores, [time_spatial_hybrid[i]/60. for i in range(len(num_cores))], marker='o', color=colors["dark green"], label="Spatial-temporal-hybrid model")
print(speed_up_spatial_hybrid)
x = np.arange(len(num_cores))
ax.plot(x, speed_up_spatial_hybrid, color=colors["purple"], marker='o', label="Speed-up")
ax.set_xticks(x)
ax.set_xticklabels(num_cores)
ax.set_ylabel("Speed-up")
ax.set_xlabel("Cores [#]")
ax2.set_xscale('log', base=2)
ax2.set_xlabel("Cores [#]")
ax2.set_xticks(num_cores)
ax2.set_xticklabels(num_cores)
ax.grid(True, which="both", ls="--", linewidth=0.5)
ax2.set_ylabel("Runtime [min]")
ax2.set_yscale("log")
fig.savefig(save_dir + "/speed_up.png", dpi=dpi)
fig2.savefig(save_dir + "/runtime.png", dpi=dpi)
handles1, labels1 = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles = handles1 + handles2
labels = labels1 + labels2
fig_leg = plt.figure(figsize=figsize)                   
fig_leg.legend(handles, labels, loc='center')      
fig_leg.tight_layout()
fig_leg.savefig(save_dir + f"/legend_speedup.png", dpi=dpi, bbox_inches='tight', transparent=True)
plt.close(fig_leg)
plt.close(fig)
