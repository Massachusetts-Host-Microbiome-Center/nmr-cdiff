#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  03 10:56:20 2022

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Synchronize multiple NMR runs using 1H spectra
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

from curveshapes import logi_c, inv_logi_c
from process import main as process

def synchronizers(runs, stretch=True, plot=True):
    """Synchronize spectra for downstream analyses.

    Parameters:
    runs -- iterable containing filepaths of each run to be synchronized; the
        first run will be the reference to which the other runs are shifted.
    stretch -- boolean flag to normalize the time scales by metabolic rate

    Returns a dictionary mapping each run base-name to a linear function that
        transforms a time vector.
    """
    if plot:
        xx = np.linspace(0, 36, num=36*10)

    # sync_functions = dict()
    for ii, filepath in enumerate(runs):
        basename = os.path.basename(filepath) # current run name

        # Process 1H stack if necessary, will be written to an excel file
        if not os.path.exists(f'{filepath}/{basename}_1H.xlsx'):
            if basename == '20210519_13CGlc':
                init = '51'
            else:
                init = '11'
            # process('1H', init, False)

        # Load "area" sheet of excel file
        areas = pd.read_excel(
            f'{filepath}/{basename}_1H.xlsx', sheet_name='area',
            engine='openpyxl', usecols=lambda x: x not in ['Scans'],
        )
        # Identify reference peak and calculate curve fit
        if "13CLeu" in basename:
            # Peak is split in 13C-Leu condition, so the shift must be adjusted
            refpk = "Isocaproate 0.76"
            tsh = 18
        else:
            refpk = "Isocaproate 0.8642"
            tsh = 10
        mask = areas['Time'] <= 36
        T = np.array(areas['Time'])[mask]
        S = np.array(areas[refpk])[mask]
        popt, pcov = scipy.optimize.curve_fit(
            logi_c, T, S, p0=[max(S), 0.3, tsh, 0]
        )
        L = popt[0]
        C = popt[3]
        popt[0] = 1
        popt[3] = 0
        c = basename
        P = inv_logi_c(np.array([0.05, 0.95]), *popt)
        if plot:
            ax = plt.axes()
            color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(xx, logi_c(xx, *popt), color=color, lw=1, label=basename)
            ax.scatter(T, (S - C)/L, marker='o', c=color,
                       facecolors=color, edgecolors='w', linewidths=0.5,
                       label=f"_{basename}")
            plt.show()
        if ii == 0:
            R = P
            yield lambda x: x
        else:
            print(f"{basename}: scaling window {P} to {R}.")
            yield lambda x: x - P[0] + R[0]
            # yield lambda x: (R[1] - R[0])/(P[1] - P[0])*(x - P[0]) + R[0]

if __name__ == "__main__":
    synchronizers([fpath for fpath in sys.argv[1:]])
