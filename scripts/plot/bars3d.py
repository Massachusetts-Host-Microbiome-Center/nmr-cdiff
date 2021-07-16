#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 10:32:45 2021

@author: apavao
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def bars3d(df, ddf=None):
    """
    Plot 3d bar plot of metabolite data (ex. SCFA, dFBA profiles)

    Parameters
    ----------
    df : pandas DataFrame
        Data to include in the plot.
        Column headers are conditions (ex. metabolite, reaction)
        Indices are timepoints
    ddf: pandas DataFrame
        Standard deviations for error bars, optional.

    Returns
    -------
    None, but displays and saves the plot.

    """
    # setup the figure and axes
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    
    t = df.index.to_numpy()
    _x = np.arange(t.size)
    m = df.columns.to_numpy()
    _y = np.arange(m.size)
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    h = df.to_numpy()
    hr = h.ravel(order='F')
    z = np.zeros_like(hr)
    w = d = 0.9
    
    if ddf is not None:
        e = ddf.to_numpy().ravel(order='F')
    
    xa = x[hr > 0]
    ya = y[hr > 0]
    za = z[hr > 0]
    ha = hr[hr > 0]
    if ddf is not None:
        ea = e[hr > 0]
    
    for i in range(_xx.shape[0]):
        ax1.plot3D(_xx[i], _yy[i], h.T[i])
        ax1.stem(_xx[i], _yy[i], h.T[i])
    #print(p)
    #p.set_sort_zpos(12)
    #print(p.do_3d_projection())
    if ddf is not None:
        for i in range(ha.size):
            ax1.plot([xa[i], xa[i]], 
                     [ya[i], ya[i]],
                     [ha[i] - ea[i], ha[i] + ea[i]],
                     '-_k', 
                     zorder=i+xa.size,
                     )
    ax1.set_title('C. diff monoassociated')
    ax1.set_xticks(_x)
    ax1.set_xticklabels(t)
    ax1.set_yticks(_y)
    ax1.set_yticklabels(m)


    
    plt.show()

def main(fn, snv, snd=None):
    """
    Main function to handle input data.

    Parameters
    ----------
    fn : string
        Filename of Excel spreadsheet.
    snv : string
        Sheet to import for values.
    snd : string
        Sheet to import for standard devs (optional).

    Returns
    -------
    None.

    """
    df = pd.read_excel(fn, sheet_name=snv, index_col=0)
    if snd is not None:
        ddf = pd.read_excel(fn, sheet_name=snd, index_col=0)
    else:
        ddf = None
    bars3d(df, ddf=ddf)
   
def valid_args(argv):
    if len(argv) < 2:
        print("Not enough arguments. Please provide excel filename and sheet name.")
        return False
    if not argv[0].endswith('.xlsx'):
        print("Please provide a valid excel spreadsheet.")
        return False
    if len(argv) > 3:
        print("Too many arguments.")
        return False
    return True
   
if __name__ == "__main__":
    args = sys.argv[1:]
    if valid_args(args):
        main(*args)