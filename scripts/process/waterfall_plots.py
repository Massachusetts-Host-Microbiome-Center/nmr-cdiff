import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import argparse
import json
import seaborn as sns

import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
import configparser
import math
from scipy.optimize import curve_fit
import trajectories
import scipy

def fit_traj(x_data, y_data, max_concentration, curve = 'logistic_C_r'):
    coeffs = trajectories.CoeffsWriter()

    # Scale and crop time series
    time = x_data
    #if tscale is not None:
    #    time = tscale(time)
    mask = (time >= 0) & (time <= 36)
    time = time[mask]
    areas = y_data[mask]

    # Perform logistic fit on products, calculate "gscale" scale factor for experiment
    signals = dict()
    curves = dict()
    curves_err = dict()

    signal = areas
    signal[trajectories.fill_mask(signal, curve)] = 0
    function = trajectories.curvedata[curve]["function"]
    p0 = trajectories.curvedata[curve]["p0"]
    scale_mask = trajectories.curvedata[curve]["scale_mask"]
    nan_mask = ~np.isnan(signal)
    popt, pcov = scipy.optimize.curve_fit(function, time[nan_mask], signal[nan_mask], p0=p0)

    amplitude = np.sum(popt[np.array(scale_mask, dtype=bool)])
    gscale = (max_concentration/amplitude)**np.array(scale_mask) # only scale L (and C)

    perr = np.sqrt(np.diagonal(pcov))
    coeffs.add_raw("", popt, perr)

    # Scale signal to concentration using gscale array
    popt *= gscale
    perr *= gscale
    coeffs.add_scaled("", popt, perr)

    return(coeffs, gscale[0])

def readPeaks(experiment_filepath):
	# read in and organize the peak position list
	peaks = pd.read_excel(experiment_filepath, header = None, sheet_name = "cfg", index_col = 0)
	peaks.columns = ['metabolite']
	peaks['colors'] = ['black'] * peaks.shape[0]
	colors = ['red', 'green', 'blue', 'turquoise', 'lightcoral', 'purple', 'orange', 'lightsteelblue', 'sienna', 'pink', 'lightgreen', 'yellow']
	for i, met in enumerate(np.unique(peaks['metabolite'])):
		peaks.loc[peaks.metabolite == met, 'colors'] = colors[i]

	return(peaks)

def readData(experiment_filepath):
	# read in and organize data
	df = pd.read_excel(experiment_filepath, sheet_name= "trace", index_col = 0)
	time = df.iloc[0, :].values
	df = df.iloc[1:]
	df.columns = time
	df[df<0] = 0
	df.index = np.array(df.index.values, dtype = float)
	return(df)

def aggPeaks(df, peaks, tolerance = 0.01):
	peak_values = peaks.index.values

	# Sum the area under the curve for each peak, where we include peaks that are within the tolerance of the target
	df_agg = pd.DataFrame()
	keep_peaks = [False] * df.shape[0]
	for peak in peak_values:
		keep = [i < (peak + tolerance) and i > (peak - tolerance) for i in df.index.values]
		df_agg[peak] = df.loc[keep, :].sum(axis = 0).values
	df_agg.index = df.columns
	df_agg = df_agg.T
	return(df_agg)

def plot(df, peaks):
	# PLOT

	ppm = df.index.values
	time = df.columns.values


	# 1. verts is a list, where every element is a list of tuples. An element is a single polygon (i.e one ppm).
	# the tuples are the (time, abundance) pairs that correspond to that ppm
	fig = plt.figure(figsize = (10, 10))
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(azim = 170, elev = 15)
	xs = df.columns.values # times, up to 48 hr

	verts = []
	print(len(verts))
	ys = df.index.values # ppms, up to about 160
	for y in ys:
		zs = np.sqrt(df.loc[y, :])
		zs[0], zs[-1] = 0, 0
		verts.append(list(zip(xs, zs)))

	poly = PolyCollection(verts, facecolors = peaks.loc[ys, 'colors'].values)
	poly.set_alpha(0.5)
	ax.add_collection3d(poly, zs=ys, zdir='y')

	ax.set_xlim(0, 100) # time 
	ax.set_ylim(0, 200) #ppms 0
	ax.set_zlim3d(0, 50)


	# manual legend
	patches = []
	for met in np.unique(peaks.metabolite.values):
		print(peaks.loc[peaks.metabolite == met, 'colors'].values[0])
		patches.append(mpatches.Patch(color=peaks.loc[peaks.metabolite == met, 'colors'].values[0], label=met))
	ax.legend(handles= patches, bbox_to_anchor = [1.1, 1.01])                  
	return(ax)

def waterfall_plot(experiment_filepath, tolerance = 0.01):
	peaks = readPeaks(experiment_filepath)
	df = readData(experiment_filepath)

	df = aggPeaks(df, peaks, tolerance)



	ax = plot(df, peaks)



def logistic(x, L, k, x0, C):
	return C + L / (1 + np.exp(-k * (x - x0)))


concentrations = {
	'Glucose': 30,
	'Acetate': 18,
	'Alanine': 15,
	'Butyrate': 7,
	'Ethanol': 2,
	'CO2': 10

}

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='waterfall_plots.py')
	parser.add_argument('--cfg_path', metavar='CONFIG', help='path to dFBA config file')
	args = parser.parse_args()
	config = configparser.ConfigParser()
	param_file = args.cfg_path
	config.read(param_file)

	project_dir = config.get('paths', 'project_dir')
	run_dir = project_dir + config.get('paths', 'run_dir')
	experiment_filepath = run_dir + config.get('paths', 'output_filename') + "_" + config.get('data', 'isotope') + ".xlsx"
	met_remove = config.get('plotting_params' ,'met_remove')
	#ax = waterfall_plot(experiment_filepath)
	#plt.savefig(run_dir + "/waterfall_plot.png")


	# plot logistic curves
	areas = pd.read_excel(experiment_filepath, sheet_name="area")
	areas = areas.melt(id_vars = ['Time'])
	keep = [i not in met_remove for i in areas['variable'].values]
	areas = areas.loc[keep, :]



	areas_fit = []
	for met in areas['variable'].unique():

				print(met)
				keep = areas['variable'] == met
				df = areas.loc[keep, ]
				keep = [not math.isnan(i) for i in df['value'].values]
				df = df.loc[keep, :]

				# Fit logistic curve to the data
				x_data  = df['Time'].values
				y_data = df['value'].values

				try:
					coeffs, gscale = fit_traj(x_data, y_data, max_concentration = concentrations[met])
					L, k, x0, C = [float(i) for i in coeffs.sopt[0][1: ]]
					#popt, pcov = curve_fit(logistic, x_data, y_data, p0=[max(y_data), 1, np.median(x_data)])
					#L, k, x0 = popt

					x_fit = np.linspace(min(x_data), max(x_data), 500)
					y_fit = logistic(x_fit, L, k, x0, C)

					tmp = pd.DataFrame({"Time": x_fit, "value": y_fit, "variable": [met] * len(x_fit)})
					areas_fit.append(tmp)

					keep = areas['variable'] == met
					print("Gscale: ", gscale)
					areas.loc[keep, 'value'] = areas.loc[keep, 'value'].values * gscale
					
				except:
					print("Could not fit: ", met)

	plt.figure(figsize=(10, 6))  # Adjust width and height in inches
	sns.scatterplot(areas, x = "Time", y = "value", hue = "variable")

	df = pd.concat(areas_fit)
	ax = sns.lineplot(df, x = "Time", y = "value", hue = "variable")
	ax.legend(bbox_to_anchor = [1.27, 1.03])
	plt.tight_layout()
	plt.savefig(run_dir + "/logistic_curves.png", 	dpi = 100)
