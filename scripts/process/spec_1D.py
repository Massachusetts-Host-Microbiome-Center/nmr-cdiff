#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 12:00:24 2021

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Process raw Bruker FID file into fourier transformed signal vector
 - Use python 3.8+ for best results
 - See /nmr-cdiff/venv/requirements.txt for dependencies
 
Copyright 2021 Massachusetts Host-Microbiome Center

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

import os, subprocess, datetime, time, glob, sys
import nmrglue
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline

SCDIR = os.path.dirname(__file__)   # location of script

def process_spec(loc, item, isotope):
    """Get process parameters using nmrglue and write NMRpipe script
    
    Arguments:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor
    
    Returns: ppm correction to calibrate the x axis
    """

    ## CONVERT and LOAD SPECTRUM ##
    subprocess.run(
        ["csh", SCDIR + "/bruker_" + isotope + ".com", item, "f"],
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    
    p0s = 0.0
    p1s = 0.0
    if isotope == '13C':
        p0s = 81.0#76.0 #91.0 # 95
        p1s = 84.0#84.0 # 68.0  #66
    elif isotope == '1H':
        p0s = 168.0#-15.0 #-90.0
        p1s = 35.0#32,0 #6.0
    elif isotope == '31P':
        p0s = 15.0
        p1s = -117.0
    
    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/ft_proc.com", item + "_f"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_f_1D.ft")
    
    #INSERT LOOP#
    p0t = p0s
    p1t = p1s
    while True:
        print(p0t, p1t)
        dic1, data1 = nmrglue.process.pipe_proc.ps(dic, data, p0=p0t, p1=p1t)
        uc = nmrglue.pipe.make_uc(dic1, data1, dim=0)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(uc.ppm_scale(), data1.real, lw=0.5)
        ax.set_xlim(200.0, 0.0)
        plt.show()
        text = input("")
        if text in ['done', 'exit', 'quit', 'q', '']:
            break
        vals = text.split()
        try:
            pa = float(vals[0])
            pb = float(vals[1])
            p0t = pa
            p1t = pb
        except ValueError:
            print("Invalid shift string.")
    #END INSERT#
    
    p0s = p0t
    p1s = p1t
    
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=p0s, p1=p1s)
    nmrglue.pipe.write(f"{loc}/{item}/{item}_f_1D.ft.ps", dic, data, overwrite=True)
    
    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/bl_proc.com", item + "_f"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## CALIBRATE PPM SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_f_1D.ft.ps.bl")
    # peaks = nmrglue.analysis.peakpick.pick(data, 1000.0)
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(uc.ppm_scale(), data.real, lw=0.5)
    ax.set_xlim(200.0, 0.0)
    print("Find the ppm shift of the peak you wish to calibrate.")
    print("Note the x-coordinate, then close the window.")
    plt.show()
    calib = input("Type the ppm shift of the peak you wish to calibrate: ")
    refpk = input("Type the reference ppm shift for the peak: ")
    cf = float(refpk) - float(calib)
    
    ## SHIFT DATA ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_f_1D.ft.ps.bl")
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    total_ppm = [i + cf for i in uc.ppm_scale()]
    
    ## CLEAN UP DATA ##
    ppm = []
    trace = []
    for i, signal in enumerate(data.real):
        ## SKIP TRIM: can trim further in MATLAB ##
        # if total_ppm[i] >= 0.0 and total_ppm[i] < 201.0:
        ppm.append(total_ppm[i])
        trace.append(signal)
    
    return ppm, trace
 
def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def subtract(loc, item, isotope):
    ## GET experimental spec
    ppm, trace = process_spec(loc, item, isotope)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ppm, trace, lw=0.5)
    ax.set_xlim(54.5, 52.5)
    plt.show()
    
    ## GET 13C-Ala spec
    refdir = "/Users/apavao/Desktop/NMR runs/May252021_C13_HRMAS_Cells"
    os.chdir(refdir)
    subprocess.run(
        ["csh", SCDIR + "/bruker_" + isotope + "_2.com", "13", "f"],
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    os.chdir(refdir + "/13")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/ft_proc_2.com", "13_f"])
    time.sleep(1.)
    os.chdir(refdir)
    dic, data = nmrglue.pipe.read(refdir + "/13/13_f_1D.ft")
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=81.0, p1=84.0)
    nmrglue.pipe.write(refdir + "/13/13_f_1D.ft.ps", dic, data, overwrite=True)
    os.chdir(refdir + "/13")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/bl_proc.com", "13_f"])
    time.sleep(1.)
    os.chdir(refdir)
    dic, data = nmrglue.pipe.read(refdir + "/13/13_f_1D.ft.ps.bl")
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    refppm = [i + 2.299 for i in uc.ppm_scale()]
    reftrace = [i for i in data.real]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(refppm, reftrace, lw=0.5)
    ax.set_xlim(54.5, 52.5)
    plt.show()
    
    refpp2 = [i + 53.643 - 53.7 for i in refppm]
    ## GET regions for normalization
    def getscale(exppm, extrc, rfppm, rftrc, lb, ub):
        exp_ppm_norm = []
        exp_trc_norm = []
        ref_ppm_norm = []
        ref_trc_norm = []
        for i, shf in enumerate(exppm):
            if shf > lb and shf < ub:
                exp_ppm_norm.append(shf)
                exp_trc_norm.append(extrc[i])
        for i, shf in enumerate(rfppm):
            if shf > lb and shf < ub:
                ref_ppm_norm.append(shf)
                ref_trc_norm.append(rftrc[i])
        exp_area = integrate.trapezoid(exp_trc_norm, exp_ppm_norm)
        ref_area = integrate.trapezoid(ref_trc_norm, ref_ppm_norm)
        return exp_area/ref_area
    
    ref_scal = getscale(ppm, trace, refppm, reftrace, 53.7, 53.75)
    xexp = []
    yexp = []
    xref = []
    yref = []
    for i, shf in enumerate(ppm):
        if shf > 50 and shf < 56:
            xexp.append(shf)
            yexp.append(trace[i])
    for i, shf in enumerate(refppm):
        if shf > 50 and shf < 56:
            xref.append(shf)
            yref.append(reftrace[i])
    for l in [xexp, yexp, xref, yref]:
        l.reverse()
    expfit = UnivariateSpline(xexp, yexp)
    reffit = UnivariateSpline(xref, yref)
    x = np.linspace(50, 56, 5000)
    y = expfit(x) - ref_scal*reffit(x)
    ref_sca2 = getscale(x, y, x, reffit(x - 53.643 + 53.7), 53.626, 53.665)
    y2 = y - 1.8*ref_sca2*reffit(x - 53.643 + 53.7)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y2, lw=0.5)
    ax.set_xlim(54.5, 52.5)
    plt.show()

def main(loc, fid, isotope='13C', verbose=True):
    """Process entire NMR stack and write to xlsx file.
    
    Arguments:
    loc -- path to folder containing traces
    ids -- list of IDs
    initial -- ID of initial acquisition, if different from first ID \
               (to find initial timepoint)
    verbose -- if True, print status messages
    """
    message("Begin NMR stack process.", verbose)

    ## calibrate and get correction factor ##
    message("Calibrate autophase and PPM shift...", verbose)
    
    subtract(loc, fid, isotope)
    
    return
    
    ppm, trace = process_spec(loc, fid, isotope)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ppm, trace, lw=0.5)
    ax.set_xlim(200, 0)
    plt.show()
    
    return
    ## write excel file ##
    message("Trace processing complete.", verbose)
    message("Writing full stack to Excel file.", verbose)
    stack_array = np.array(stack)
    df = pd.DataFrame(stack_array[1:], columns=stack_array[0], index=header[1])
    writer = pd.ExcelWriter(loc + '/StackSpec_' + isotope + '.xlsx', engine='xlsxwriter')
    df = df.T
    df.to_excel(writer, sheet_name='trace', startrow=1)
    wb = writer.book
    ws = writer.sheets['trace']
    ws.write(0, 0, 'trace')
    ws.write(1, 0, 'ppm')
    for i, item in enumerate(header[0]):
        ws.write(0, i + 1, item)
    if os.path.exists(loc + '/cfg.txt'):
        # Write reference ppm config
        message("Writing reference peak config.", verbose)
        cfgsheet = wb.add_worksheet('cfg')
        with open(loc + '/cfg.txt', 'r') as rf:
            for i, line in enumerate(rf):
                l = line.strip('\n').split('\t')
                cfgsheet.write(i, 0, l[0])
                cfgsheet.write(i, 1, l[1])
    else:
        message("Config file not found.", verbose)
    writer.save()
    message(f"Completed successfully. Output: {loc}/StackSpec_{isotope}.xlsx", verbose)
    
def handle_args(args):    
    if len(args) < 1:
        iso = '13C'
    elif args[0] in ['-h', '--help']:
        print("Spectral processing script.")
        print("===========================")
        print("Command line args:")
        print("\t1. Isotope. Supported: 1H, 13C (default), 31P.")
        print("\t2. Spectrum ID. Default=11.")
        return '13C', '11'
    elif args[0] in ['13C', '1H', '31P']:
        iso = args[0]
    else:
        try:
            fid = str(int(args[0]))
        except ValueError:
            print("Invalid arguments.")
            return None
    if len(args) > 1:
        try:
            fid = str(int(args[1]))
        except ValueError:
            print("Invalid arguments.")
            return None
    else:
        fid = '11'
    if len(args) > 2:
        print("Ignoring excess arguments: " + " ".join(args[2:]))
    loc = os.getcwd()
    return loc, fid, iso

if __name__ == "__main__":   
    main(*handle_args(sys.argv[1:]))