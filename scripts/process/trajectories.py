#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  22 11:35:45 2022

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Estimate concentration trajectories from HRMAS NMR time series and standard
   solutions
 - Use python 3.8+ for best results
 - See /nmr-cdiff/venv/requirements.txt for dependencies

Copyright 2022 Massachusetts Host-Microbiome Center

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from process import main as process
from standard_curves import compute_standard_curves
from curveshapes import logi, logi_c, logi_flex
from get_color import get_cmap

SCDIR = os.path.dirname(__file__)   # location of script
STANDARDS_FILE = SCDIR + '/../../data/test/20220421_13CGlc_Standards'

def find_substrate(basename):
    """Guess name, concentration of substrate from filename."""
    substrate_map = {
        # 13C-Compound: (Name, Concentration in MMM (mM)),
        '13CGlc': ('Glucose', 27.78),
        '13CPro': ('Proline', 6.96),
        '13CLeu': ('Leucine', 7.63),
        '13CAc': ('Acetate13C', 13.445)
    }
    for k, v in substrate_map.items():
        if k in basename:
            return v
    return None

def fit_trajectories(filepath, init='11', plot=False, tscale=None,
                     spath=STANDARDS_FILE):
    """Estimate concentration trajectories for HRMAS NMR time series.

    Parameters:
    filepath -- path to the directory containing the Bruker run. If an Excel
        spreadsheet with the output from process.py does not exist, the NMR
        experiment will be processed.
    init -- initial spectrum for processing script, if it will be run
    plot -- True to show a plot of the trajectories
    tscale -- a linear function mapping run time to scaled time. Recommended if
        multiple runs are being analyzed, can be calculated using
        synchronizers.py
    spath -- if the substrate requires a set of NMR standards to estimate
        concentration, the filepath to an experiment measuring those standards

    Returns a dictionary of scaled logistic coefficients for each metabolite,
    and dictionaries of lower and upper bounds for those coefficients from
    scipy.
    """

    data = {'popt': [], 'perr': [], 'sopt': [], 'serr': []}

    basename = os.path.basename(filepath) # current run name

    print("=========")
    print(f"Computing logistic trajectories for run {basename}.")
    print("=========")

    # Process stack, will be written to an excel file
    if not os.path.exists(f'{filepath}/{basename}_13C.xlsx'):
        process('13C', init)

    # Load "area" sheet of excel file
    areas = pd.read_excel(
        f'{filepath}/{basename}_13C.xlsx', sheet_name='area',
        engine='openpyxl', usecols=lambda x: x not in ['Scans'],
    )

    substrate, input_cx = find_substrate(basename)
    sareas = areas[substrate]
    T = areas['Time']
    if tscale is not None:
        T = tscale(T)

    # Substrate and products are scaled by standard curves
    if substrate == "Leucine":
        sfacts = {
            "Isocaproate": 0.633,
            "Isovalerate": 0.367,
            "Isobutyrate": 0.170,
            "2-methylbutyrate": 0.3,
        }
    elif substrate == "Proline":
        sfacts = {
            "5-aminovalerate": 1.,
        }
    elif substrate == "Acetate13C":
        sfacts = {
            "AcetateNA": 0.5,
        }
    elif substrate == "Glucose":
        print("Computing standard curves.")
        print("---------")
        sfacts = compute_standard_curves(filepath=spath)
        print("---------")
    pareas = areas[[c for c in areas.columns if c in sfacts]]
    cest = pareas.apply(lambda col: col*sfacts[col.name], axis=0)

    cest = cest.fillna(0)

    # Scale spectra
    curves = dict()
    curves_err = dict()

    if plot:
        fig = plt.figure(figsize=(2.5, 2), constrained_layout=True)
        ax = plt.axes()
        xx = np.linspace(0, 36, num=36*10)
        cmap = get_cmap()
        plt.axhline(ls=':', color=cmap['gray'], linewidth=1, zorder=1)
    # Get substrate curve
    mask = (T >= 0) & (T <= 36)
    A = np.array(sareas)[mask]
    if substrate == "Acetate13C":
        popt, pcov = scipy.optimize.curve_fit(
            logi, np.array(T[mask]), A, p0=[500, 0.2, 12]
        )
        gscale = (input_cx/popt[0])**np.array([1, 0, 0]) # only scale L and C
    else:
        popt, pcov = scipy.optimize.curve_fit(
            logi_c, np.array(T[mask]), A, p0=[1000, -0.2, 12, 400]
        )
        gscale = (input_cx/(popt[0] + popt[3]))**np.array([1, 0, 0, 1]) # only scale L and C
    perr = np.sqrt(np.diagonal(pcov))
    data['popt'].append([substrate] + [f"{ele:.3f}" for ele in popt])
    data['perr'].append([substrate] + [f"{ele:.3f}" for ele in perr])

    # Scale signal to concentration using gscale array
    popt *= gscale
    perr *= gscale
    data['sopt'].append([substrate] + [f"{ele:.3f}" for ele in popt])
    data['serr'].append([substrate] + [f"{ele:.3f}" for ele in perr])
    curves[substrate] = popt
    curves_err[substrate] = perr
    print("Logistic coefficients")
    print("---------")
    print(substrate)
    print(f" - L : {popt[0]}")
    print(f" - k : {popt[1]}")
    print(f" - x0: {popt[2]}")
    if substrate != "Acetate13C":
        print(f" - C : {popt[3]}")
        gscale = gscale[:-1] # products do not have C coefficient

    if plot:
        color = cmap[substrate]
        ax.plot(
            xx, logi_flex(xx, *popt), color=color, lw=1, label=substrate,
            zorder=2,
        )
        im = ax.scatter(
            T[mask], A*gscale[0], marker='o', color=color, facecolors=color,
            edgecolors='w', linewidths=0.5, sizes=[22 for ele in A], zorder=3,
            label=f"_{substrate}",
        )


    # Get products curves
    for i, c in enumerate(cest.columns):
        mask = (~np.isnan(cest[c]))
        A = np.array(cest[c])
        try:
            popt, pcov = scipy.optimize.curve_fit(
                logi, np.array(T[mask]), A, p0=[500, 0.2, 12]
            )
        except:
            print(f"Optimal parameters not found for compound {c} and run {basename}.")
            continue
        perr = np.sqrt(np.diagonal(pcov))
        data['popt'].append([c] + [f"{ele:.3f}" for ele in popt])
        data['perr'].append([c] + [f"{ele:.3f}" for ele in perr])
        if c in ["Isocaproate", "Isovalerate", "5-aminovalerate"]:
            gscale = (input_cx*sfacts[c]/popt[0])**np.array([1, 0, 0])
        popt *= gscale
        perr *= gscale
        data['sopt'].append([c] + [f"{ele:.3f}" for ele in popt])
        data['serr'].append([c] + [f"{ele:.3f}" for ele in perr])
        curves[c] = popt
        curves_err[c] = perr

        print(c)
        print(f" - L : {popt[0]}")
        print(f" - k : {popt[1]}")
        print(f" - x0: {popt[2]}")

        ## Add curves for other oxidative BCAA fermentations
        if c == "Isovalerate":
            for subst, prod in [("Valine", "Isobutyrate"), ("Isoleucine", "2-methylbutyrate")]:
                iscale = (sfacts[prod]/sfacts[c])**np.array([1, 0, 0])
                popt_i = popt*iscale
                perr_i = perr*iscale
                data['sopt'].append([prod] + [f"{e:.3f}" for e in popt_i])
                data['serr'].append([prod] + [f"{e:.3f}" for e in perr_i])
                curves[prod] = popt_i
                curves_err[prod] = perr_i
                vscale = (-1)**np.array([0, 1, 0])
                popt_v = popt_i*vscale
                data['sopt'].append([subst] + [f"{e:.3f}" for e in popt_v])
                data['serr'].append([subst] + [f"{e:.3f}" for e in perr_i])
                curves[subst] = popt_v
                curves_err[subst] = perr_i

        if plot:
            color = cmap[c]
            ax.plot(
                xx, logi_flex(xx, *popt), color=color, lw=1, label=c, zorder=2,
            )
            im = ax.scatter(
                T[mask], A*gscale[0], marker='o', color=color, facecolors=color,
                edgecolors='w', linewidths=0.5, label=f"_{c}", zorder=3,
                sizes=[22 for ele in A],
            )
    if plot:
        afont = {'fontname': 'Arial', 'size': 7}
        ax.set_xlabel("Time (h)", **afont)
        ax.set_ylabel("Estimated Concentration (mM)", **afont)
        ax.set_xlim((0, 36))
        if substrate == "Glucose":
            ax.set_yticks([0, 10, 20, 30])
            ax.set_yticks(list(range(-2, 31)), minor=True)
        elif substrate == "Proline":
            ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
            ax.set_yticks(list(np.arange(-0.5, 8, 0.5)), minor=True)
        elif substrate == "Leucine":
            ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
            ax.set_yticks(list(np.arange(-0.5, 9, 0.5)), minor=True)
        ax.set_xticks([0, 12, 24, 36])
        ax.set_xticks(list(range(37)), minor=True)
        ax.set_xticklabels(ax.get_xticks(), **afont)
        ax.set_yticklabels(ax.get_yticks(), **afont)
        ax.xaxis.set_tick_params(width=0.5)
        ax.yaxis.set_tick_params(width=0.5)
        plt.setp(ax.spines.values(), linewidth=0.5)
        # plt.legend()
        # plt.show()
        plt.savefig(f"{filepath}/{basename}_pan.svg")

    for l in ['p', 's']:
        with open(f"{filepath}/{basename}_logistic{l}.csv", "w") as wf:
            opt = data[f"{l}opt"]
            err = data[f"{l}err"]
            wf.write("Met,L,k,x0,C\n")
            for i, (o, e) in enumerate(zip(opt, err)):
                segs = [o[0]] + [f"{k}pm{m}" for (k, m) in zip(o[1:], e[1:])]
                wf.write(','.join(segs) + ','*(5-len(segs)) + '\n')

    return curves, curves_err #curves_lb, curves_ub

def call_fit(args):
    n_args = len(args)
    if n_args > 1:
        print("Too many arguments. Please give only the run directory, or no "
              + "arguments to use the current working directory.")
    else:
        if n_args == 0:
            filepath = os.getcwd() # path to working directory
        else:
            filepath = args[0]
        return fit_trajectories(filepath, plot=True, scaled=True)

if __name__ == "__main__":
    call_fit(sys.argv[1:])
