#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 12:00:24 2021

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Process raw Bruker NMR files into Excel spreadsheet for downstream analysis
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

SCDIR = os.path.dirname(__file__)   # location of script

def calibrate_process(loc, item, acq1, isotope):
    """Get process parameters using nmrglue and write NMRpipe script
    
    Arguments:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor
    
    Returns: ppm correction to calibrate the x axis
    """
    ## GET INITIAL TIMESTAMP ##
    time0 = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{acq1}/pulseprogram"))
    
    ## CONVERT and LOAD SPECTRUM ##
    subprocess.run(
        ["csh", SCDIR + "/bruker_" + isotope + "_2.com", item, "calibr"],
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
    subprocess.run(["csh", SCDIR + "/ft_proc.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft")
    
    #INSERT LOOP#
    # p0t = p0s
    # p1t = p1s
    # for i in range(8):
    #     for j in range(8):
    #         p0t = 75.0 + 2.0*i
    #         p1t = 70.0 + 2.0*j
    #         print(p0t, p1t)
    #         dic1, data1 = nmrglue.process.pipe_proc.ps(dic, data, p0=p0t, p1=p1t)
    #         uc = nmrglue.pipe.make_uc(dic1, data1, dim=0)
    #         fig = plt.figure()
    #         ax = fig.add_subplot(111)
    #         ax.plot(uc.ppm_scale(), data1.real, lw=0.5)
    #         plt.show()
    #END INSERT#
    
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=p0s, p1=p1s)
    nmrglue.pipe.write(f"{loc}/{item}/{item}_calibr_1D.ft.ps", dic, data, overwrite=True)
    
    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/bl_proc.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## CALIBRATE PPM SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft.ps.bl")
    # peaks = nmrglue.analysis.peakpick.pick(data, 1000.0)
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(uc.ppm_scale(), data.real, lw=0.5)
    print("Find the ppm shift of the peak you wish to calibrate.")
    print("Note the x-coordinate, then close the window.")
    plt.show()
    calib = input("Type the ppm shift of the peak you wish to calibrate: ")
    refpk = input("Type the reference ppm shift for the peak: ")
    cf = float(refpk) - float(calib)
    
    return cf, time0, [p0s, p1s]

def process_trace(loc, item, cf, phases, isotope):
    """Process a single 1D 13C-NMR trace.
    
    Arguments:
    loc -- path to folder containing traces
    item -- ID of trace to process
    cf -- correction factor for ppm shift
    
    Returns:
    ppm --  list of ppm values
    trace -- list of signal intensities
    tim --  timestamp of trace, calculated from acquisition data
    phases -- P0, P1 phase shift
    """
    ## GET TIMESTAMP ##
    time0 = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{item}/pulseprogram"))
    time1 = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{item}/acqus"))
    tim = time0 + (time1 - time0)/2
    
    ## CONVERT BRUK to PIPE ##
    subprocess.run(["csh", SCDIR + "/bruker_" + isotope + "_2.com", item],        
                   stdout=subprocess.DEVNULL, 
                   stderr=subprocess.DEVNULL)
    
    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/ft_proc.com", item])
    time.sleep(1.)
    os.chdir(loc)
    
    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_1D.ft")
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=phases[0], p1=phases[1])
    nmrglue.pipe.write(f"{loc}/{item}/{item}_1D.ft.ps", dic, data, overwrite=True)
    
    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/bl_proc.com", item])
    time.sleep(1.)
    os.chdir(loc)
    
    ## SHIFT DATA ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_1D.ft.ps.bl")
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
    
    return ppm, trace, tim
 
def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, isotope='13C', verbose=True):
    """Process entire NMR stack and write to xlsx file.
    
    Arguments:
    loc -- path to folder containing traces
    ids -- list of IDs
    initial -- ID of initial acquisition, if different from first ID \
               (to find initial timepoint)
    verbose -- if True, print status messages
    """
    message("Begin NMR stack process.", verbose)
    ids = [str(item) for item in ids]
    if initial == None:
        initial = ids[0]
    
    header = [[], []]
    stack = []
    
    ## calibrate and get correction factor ##
    message("Calibrate autophase and PPM shift...", verbose)
    cf, time0, phases = calibrate_process(loc, ids[0], initial, isotope)
    
    ## process all traces ##
    message("Begin processing " + isotope + " traces.", verbose)
    for i, item in enumerate(ids):
        message(f"Processing file {item}...", verbose)
        ppm, trace, timestamp = process_trace(loc, item, cf, phases, isotope)
        delta_t = (timestamp - time0).total_seconds()/3600.
        if i == 0:
            stack.append(ppm)
        if ppm != stack[0]:
            print("PPM mismatch. Abort.")
            return
        header[0].append(item)
        header[1].append(delta_t)
        stack.append(trace)
    
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

def detect_spectra(isotope='13C', init=11):
    """Returns a list of spectra IDs for the given isotope."""
    iso_spectra = []
    dirs = glob.glob('*/')
    st = [str(y) for y in sorted([int(x[:-1]) for x in dirs])]
    for fn in st:
        if int(fn) >= int(init):
            try:
                dic, data = nmrglue.bruker.read(fn)
            except OSError:
                break
            udic = nmrglue.bruker.guess_udic(dic, data)
            if udic[0]['label'] == isotope and udic[0]['encoding'] == 'direct':
                iso_spectra.append(fn)
    return iso_spectra

def main(iso, init):
    """Call stack processing function."""
    print("Welcome to NMR processing pipeline.\n")
    if iso not in ("1H", "13C", "31P"):
        print("Isotope not recognized.")
        return
    loc = os.getcwd()
    indices = detect_spectra(isotope=iso, init=init)
    print("Beginning calibration with following parameters:")
    print("\t- Directory: " + loc)
    print(f"\t- {iso} spectrum IDs: " + ', '.join(indices))
    print("\t- Initial timepoint reference: " + init)
    process_stack(loc, indices, initial=init, isotope=iso, verbose=True)
    
def handle_args(args):    
    if len(args) < 1:
        iso = '13C'
    elif args[0] in ['-h', '--help']:
        print("Spectral processing script.")
        print("===========================")
        print("Command line args:")
        print("\t1. Isotope. Supported: 1H, 13C (default), 31P.")
        print("\t2. Initial spectrum. Timestamp used as t0. Default=11.")
        return '13C', '11'
    elif args[0] in ['13C', '1H', '31P']:
        iso = args[0]
    else:
        try:
            init = str(int(args[0]))
        except ValueError:
            print("Invalid arguments.")
            return None
    if len(args) > 1:
        try:
            init = str(int(args[1]))
        except ValueError:
            print("Invalid arguments.")
            return None
    else:
        init = '11'
    if len(args) > 2:
        print("Ignoring excess arguments: " + " ".join(args[2:]))
    return iso, init

if __name__ == "__main__":   
    main(*handle_args(sys.argv[1:]))