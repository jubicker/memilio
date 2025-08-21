import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('default')
fontsize = 10
matplotlib.rcParams.update({
	'font.size': fontsize,
'axes.titlesize': fontsize * 1.2,
 	'axes.labelsize': fontsize,
 	'xtick.labelsize': fontsize * 0.8,
 	'ytick.labelsize': fontsize * 0.8,
	'legend.fontsize': fontsize * 0.8,
	'font.family': 'DejaVu Sans'
})
plt.rcParams['lines.linewidth'] = 1

colors = {"S": "#3868b0", "I": "#a51919", "R": "#7B7B7B"}
other_colors = ['orange', 'blue', 'red', 'green', 'gray', 'black', 'teal', 'brown']
dpi = 300
sim_results = '/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/'
save_folder = '/home/bick_ju/Documents/TemporalHybridModel/visualizations/'
