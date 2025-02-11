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






def plotLogisticCurve(nmr_run_folder, metabolite=["AlanineN15", "Alanine"], run_names=[], pvals=None, title=None):
    fig = plt.figure(figsize=(8, 6))
    marker_styles = {name: marker for name, marker in zip(run_names, ['s', '^', 'o', 'D', 'x', 'P', '*', 'H', 'v', '<'])}
    color_map = {
        metabolite[0]: 'blue',
        metabolite[1]: 'orange',
        # Add more metabolites and their corresponding colors as needed
    }

    log_params = {}
    legend_labels = set()  # Set to track added legend labels
    data_max = 0
    data_min = 0
    for run_name in run_names:
        folder = f"{nmr_run_folder}/{run_name}/results/"
        df = pd.read_csv(folder + "time_norm_areas.csv")
        log_params[run_name] = {}

        for met in metabolite:
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
            log_params[run_name][met] = popt
            print(popt)
            
            x_fit = np.linspace(time.min(), time.max(), 100)
            y_fit = logistic_func(x_fit, *popt)
            
            marker = marker_styles.get(run_name, 'o')  # Default to circle if run_name not found
            color = color_map.get(met, 'gray')  # Default to gray if metabolite not found
            
            # Plot the scatter and fit
            plt.scatter(time, data, marker=marker, color=color)
            plt.plot(x_fit, y_fit, color=color)

            if data.max() > data_max:
                data_max = data.max()
            if data.min() < data_min:
                data_min = data.min()

            # Add to legend only if not already added
            if met not in legend_labels:
                plt.plot([], [], label=met, color=color)  # Dummy plot for legend
                legend_labels.add(met)
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
        plt.annotate("Rank sum test", (5, y_start ), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color="black")
        for i, param in enumerate(pvals.index):
            p_value = pvals.loc[param, 'pvalue']  
            plt.annotate(f"{param}: p={p_value:.3f}", (5, y_start -1/5 - i/5), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color="black")

    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Concentrations (projected to standard curve)")
    plt.legend()
    
    return pd.DataFrame(log_params), fig





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nmr_run_folder', type = str)
    parser.add_argument('--output_folder', type = str)
    parser.add_argument('--run_names','--list', nargs='+', required=True)
    parser.add_argument('--plot_title', type = str)
    parser.add_argument('--metabolites_to_compare',  nargs='+', required=True)
    args = parser.parse_args()


    # 1. Get original plot w/ logistic curve parameters
    run_names = args.run_names
    print(run_names)
    log_params, fig = plotLogisticCurve(nmr_run_folder = args.nmr_run_folder, run_names = run_names, metabolite = args.metabolites_to_compare)
    print(log_params) # in order, K (slope), x0 (midpoint), L (asymptote)


    # 2. parse logistic curve parameters:
    long_form_log_params = []
    for run_name, met_dict in log_params.items():
        for met, params in met_dict.items():
            long_form_log_params.append({
                'Run Name': run_name,
                'Metabolite': met,
                'K (slope)': params[0],
                'x0 (midpoint)': params[1],
                'L (asymptote)': params[2]
            })
    long_form_log_params_df = pd.DataFrame(long_form_log_params)

    # 3. Rank sum test on logistic curve parameters
    group1_mask = long_form_log_params_df["Metabolite"] == args.metabolites_to_compare[0]
    group2_mask = long_form_log_params_df["Metabolite"] == args.metabolites_to_compare[1]

    pvals = []
    params = ["K (slope)", "x0 (midpoint)", "L (asymptote)"]
    for param in params:
        group1 = long_form_log_params_df.loc[group1_mask, param]
        group2 = long_form_log_params_df.loc[group2_mask, param]
        s, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        print(f"{param}: p-value = {p_value}")
        pvals.append(p_value)

    pvals = pd.DataFrame(pvals, index = params, columns = ["pvalue"])
    pvals

    # 4. Plot the trajectories with the p-values
    log_params, fig = plotLogisticCurve(nmr_run_folder = args.nmr_run_folder, metabolite=args.metabolites_to_compare, run_names = run_names, pvals = pvals, title = args.plot_title)
    fig.savefig(f"{args.output_folder}/trajectories_with_pvalues.png")