import argparse
import datetime
import json
import os

import cobra as cb
from matplotlib import pyplot as plt
import matplotlib.colors
import numpy as np
import openpyxl as xl
import pandas as pd

from curveshapes import Metabolite
from standard_curves import compute_standard_curves
from synchronize import synchronizers
from trajectories import fit_trajectories
from get_color import get_cmap
import dfba
import seaborn as sns
from importlib import reload
reload(dfba)
import curveshapes
reload(curveshapes)
import trajectories

from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu


def logistic_func(x, a, b, c):
    return c / (1 + np.exp(-a * (x - b)))



def plotLogisticCurve(nmr_run_folder, met="Alanine", run_names=[], run_metadata_shape=[], run_metadata_color=[], pvals=None, title=None):
    fig = plt.figure(figsize=(8, 6))
    marker_styles = {name: marker for name, marker in zip(run_metadata_shape, ['s', '^', 'o', 'D', 'x', 'P', '*', 'H', 'v', '<'])}
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']  # List of colors
    color_map = {run_metadata_color[i]: colors[i] for i in range(min(len(run_metadata_color), len(colors)))}

    log_params = {}
    legend_labels = set()  # Set to track added legend labels
    data_max = 0
    data_min = 0
    for i, run_name in enumerate(run_names):
        folder = f"{nmr_run_folder}/{run_name}/results/"
        df = pd.read_csv(folder + "time_norm_areas.csv")


        data = df[met].values
        time = df['Time'].values
        
        # Create a mask for non-NaN values
        mask = ~np.isnan(data) & ~np.isnan(time)
        data = data[mask]
        time = time[mask]

        if len(data) == 0 or len(time) == 0:
            print(f"No valid data for {met} in {run_name}. Skipping...")
            continue  # Skip this metabolite if there's no valid data

        popt, pcov = curve_fit(logistic_func, time, data)
        log_params[run_name] = popt
        print(popt)
        
        x_fit = np.linspace(time.min(), time.max(), 100)
        y_fit = logistic_func(x_fit, *popt)
        
        marker = marker_styles.get(run_metadata_shape[i], 'o')  # Default to circle if run_name not found
        color = color_map.get(run_metadata_color[i], 'gray')  # Default to gray if metabolite not found
        
        # Plot the scatter and fit
        plt.scatter(time, data, marker=marker, color=color)
        plt.plot(x_fit, y_fit, color=color)

        if data.max() > data_max:
            data_max = data.max()
        if data.min() < data_min:
            data_min = data.min()

        # Add to legend only if not already added
        
        if run_metadata_color[i] not in legend_labels:
            plt.plot([], [], label=run_metadata_color[i], color=color)  # Dummy plot for legend
            legend_labels.add(run_metadata_color[i])
    # Add marker legend
    for name, marker in marker_styles.items():
        plt.plot([], [], marker=marker, color='black', label=name)  # Dummy plot for marker legend

    # Annotate p-values if available
    if pvals is not None:
        y_min, y_max = data_min, data_max
        y_range = y_max - y_min
        print("Ymax: ", y_max)
        y_start = y_max - (y_range / 3)  # Start at half the existing range
        print("Ystart: ", y_start)
        plt.annotate("Rank sum test", (2, y_start ), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color="black")
        for i, param in enumerate(pvals.index):
            p_value = pvals.loc[param, 'pvalue']  
            plt.annotate(f"{param}: p={p_value:.3f}", (2, y_start -1/2 - i/2), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color="black")

    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Concentrations (projected to standard curve)")
    plt.legend()
    log_params = pd.DataFrame(log_params)
    log_params.index = ["K", "x0", "L"]
    return pd.DataFrame(log_params), fig



def get_pvalues(log_params, group1_mask, group2_mask):
    # pvalue
    params = ["K", "x0", "L"]
    pvals = []
    for param in params:
        group1 = log_params.loc[param, group1_mask]
        group2 = log_params.loc[param, group2_mask]
        s, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        print(f"{param}: p-value = {p_value}")
        pvals.append(p_value)
    pvals = pd.DataFrame(pvals, index = params, columns = ["pvalue"])
    return(pvals)






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nmr_run_folder', type = str)
    parser.add_argument('--output_folder', help='path to dFBA config file')
    parser.add_argument('--datasets', nargs='+', required=True)
    parser.add_argument("--mets", nargs='+', required=True)

    args = parser.parse_args()

    
    # 1. Get original plot w/ logistic curve parameters
    run_names = ["cdiff_glucose_13C_mgh_122420", "cdiff_glucose_13C_mgh_051921", "cdiff_glucose_13C_mgh_101322", "cdiff_mannitol_13C_mgh_030323", "cdiff_glucose_13C_umass_100424", 
                "cdiff_mannitol_13C_umass_110824", "cdiff_mannitol_13C_umass_112224"]


    # TODO: adjust to make more general with input arguments
    print(run_names)
    run_substrate = np.array([None] * len(run_names))
    run_substrate[["glucose" in i for i in run_names]] = "Glucose"
    run_substrate[["mannitol" in i for i in run_names]] = "Mannitol"

    run_institute = np.array([None] * len(run_names))
    run_institute[["umass" in i for i in run_names]] = "UMass"
    run_institute[["mgh" in i for i in run_names]] = "MGH"
    group1_mask = run_substrate == "Glucose"
    group2_mask = run_substrate == "Mannitol"

    for met in args.mets:
        log_params, fig = plotLogisticCurve(nmr_run_folder = args.nmr_run_folder, run_names = run_names, met = met, run_metadata_color=run_substrate, run_metadata_shape=run_institute, title = met)
        pvals = get_pvalues(log_params, group1_mask, group2_mask)
        log_params, fig = plotLogisticCurve(nmr_run_folder = args.nmr_run_folder, run_names = run_names, met = met, run_metadata_color=run_substrate, run_metadata_shape=run_institute, title = met, pvals = pvals)
        fig.savefig(f"{args.output_folder}/trajectories_with_pvalues_{met}.png")

