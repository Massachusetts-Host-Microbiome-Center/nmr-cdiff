#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:01:56 2021

Read NMR signal data from excel file

Give filename, x column name, y column name in command line.

@author: apavao
"""

import sys
import matplotlib.pyplot as plt
import pandas as pd

def main(fn, xname, yname):
    try:
        df = pd.read_excel(fn, sheet_name='trace', header=0, index_col=0, skiprows=1)
    except ImportError:
        print("Error reading excel file. Please install xlrd.")
        return
    if xname == 'index' or xname == 'ppm':
        xdata = df.index
    else:
        xdata = df[xname]
    print([x for x in df.columns])
    ydata = df[yname]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xdata, ydata, lw=0.5)
    ax.set_xlim(200, 0) # 53.8, 53.0
    plt.show()

if __name__ == '__main__':
    args = sys.argv[1:]
    main(*args)
