import numpy as np
import matplotlib.pyplot as plt
from settings import *

# Data
x = np.linspace(0, 5, 50)

# curves = [
#     0.9 * x**2 + 0.2 * x,
#     0.8 * x**2 + 0.15 * x,
#     0.5 * x**2 + 0.1 * x,
#     0.35 * x**2 + 0.05 * x,
#     0.25 * x**2,
#     0.18 * x**2 + 0.02,
#     0.12 * x**2 + 0.03
# ]

curves = [
    0.9 * x**2 + 0.2 * x,
    0.5 * x**2 + 0.1 * x,
    0.25 * x**2]

# colors = [colors['dark teal'], colors['purple'], colors['middle blue'], colors['dark blue'], colors["very dark blue"], colors['teal'], colors['middle green']]

# Figure
fig, ax = plt.subplots(figsize=(2.5,1.7))

ax.plot(x, curves[1], color="#9C0202", linewidth=1)
ax.fill_between(x, curves[0], curves[2], color="#9C0202", alpha=0.4, lw=0)

# # Plot curves
# for y, c in zip(curves, colors):
#     noise = np.random.normal(0, 0.0, size=len(x)) * (x / max(x))  # small growing noise
#     ax.plot(x, y + noise, color=c, linewidth=0.9)
    

ax.set_xlim(-0.2, 15)
ax.set_ylim(-0.2, 45)

ax.axis("off")
# Arrow axes
ax.annotate("", xy=(15,0), xytext=(-0.2,0),
            arrowprops=dict(arrowstyle="-|>", lw=1, color="black", mutation_scale=8))

ax.annotate("", xy=(0,45), xytext=(0,-1),
            arrowprops=dict(arrowstyle="-|>", lw=1, color="black", mutation_scale=8))

fig.savefig("example_curves.png", dpi=dpi, bbox_inches='tight')

plt.close(fig)
