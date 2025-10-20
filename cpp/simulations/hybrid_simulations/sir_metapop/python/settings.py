import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import numpy as np
import cycler

fontsize = 8
plt.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize * 1,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize * 0.8,
    'ytick.labelsize': fontsize * 0.8,
    'legend.fontsize': fontsize * 0.8,
    'font.family': "Arial"
})

plt.style.use('default')

dpi = 300

colors = {"dark blue": "#2A2F8D", "middle blue": "#418EC1", "light blue": "#8FCAF1", "teal": "#5AA18E", "light teal": "#AFE5DB", "dark green": "#3A6713", "middle green": "#8AA72B", "light green": "#BAE772", "purple": "#7F4390", "rose": "#D95E77", "red": "#C62A2A", "dark red": "#872C12", "brown": "#A3682A", "orange": "#EF8D25", "yellow": "#F9CD20", "light yellow": "#FBE960", "black": "#000000", "dark grey": "#656565", "middle grey": "#A4A4A4", "light grey": "#D7D7D7", "white": "#FFFFFF"}

compartment_colors = {
    "S": colors["dark grey"],
    "I": colors["dark red"],
    "R": colors["middle green"]
}
