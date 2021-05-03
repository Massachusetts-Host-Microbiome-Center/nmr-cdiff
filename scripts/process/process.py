#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 12:00:24 2021

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Process raw Bruker NMR files into Excel spreadsheet for downstream analysis
 - Use python 3.8+ for best results
 - See /nmr-cdiff/venv/requirements.txt for dependencies
"""

import os, subprocess, datetime, time
import nmrglue
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCDIR = os.path.dirname(__file__)   # location of script

def find_extreme_pks(locs):
    # find lower peak
    locs = locs[::-1]
    subpks_l = [locs[0]]
    for i in range(1, len(locs)):
        if locs[i] - locs[i - 1] < 0.5:
            subpks_l.append(locs[i])
        else:
            break
    subpks_h = [locs[-1]]
    for i in range(1, len(locs)):
        if locs[-1*i] - locs[-1*i - 1] < 0.5:
            subpks_h.append(locs[-1*i - 1])
        else:
            break
    return [sum(subpks_l)/len(subpks_l), sum(subpks_h)/len(subpks_h)]

def calibrate_process2(loc, item, acq1):
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
        ["csh", SCDIR + "/fid_bruk_13C_calib.com", item],
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL)
    
    p0s = 91.0
    p1s = 68.0
    
    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/fid_aph_13C_ft.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft")
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=p0s, p1=p1s)
    nmrglue.pipe.write(f"{loc}/{item}/{item}_calibr_1D.ft.ps", dic, data, overwrite=True)
    
    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/fid_aph_13C_bl.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    
    ## CALIBRATE PPM SHIFT ##
    ref = []
    if os.path.exists(loc + '/cfg.txt'):
        with open(loc + '/cfg.txt', 'r') as rf:
            for line in rf:
                l = line.strip('\n').split('\t')
                ref.append(float(l[0]))
    else:
        ref_str = input("Please enter reference ppm for calibration peak: ")
        ref.append(float(ref_str))
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft.ps.bl")
    # peaks = nmrglue.analysis.peakpick.pick(data, 1000.0)
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(uc.ppm_scale(), data.real, lw=0.5)
    for pk in ref:
        ax.axvline(x=pk, label=pk, ls=':')
    # plt.legend()
    plt.show()
    calib = input("Type the ppm shift of the peak you wish to calibrate: ")
    refpk = input("Type the reference ppm shift for the peak: ")
    plt.close('all')
    cf = float(refpk) - float(calib)
    
    return cf, time0, [p0s, p1s]

def process_trace2(loc, item, cf, phases):
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
    subprocess.run(["csh", SCDIR + "/fid_bruk_13C.com", item],        
                   stdout=subprocess.DEVNULL, 
                   stderr=subprocess.DEVNULL)
    
    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/fid_aph_13C_ft.com", item])
    time.sleep(1.)
    os.chdir(loc)
    
    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_1D.ft")
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=phases[0], p1=phases[1])
    nmrglue.pipe.write(f"{loc}/{item}/{item}_1D.ft.ps", dic, data, overwrite=True)
    
    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/fid_aph_13C_bl.com", item])
    time.sleep(1.)
    os.chdir(loc)
    
    ## SHIFT DATA ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_1D.ft.ps.bl")
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    total_ppm = [i + cf for i in uc.ppm_scale()]
    
    ## TRIM DATA ##
    ppm = []
    trace = []
    for i, signal in enumerate(data.real):
        if total_ppm[i] >= 0.0 and total_ppm[i] < 201.0:
            ppm.append(total_ppm[i])
            trace.append(signal)
    
    return ppm, trace, tim
 
def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, verbose=True):
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
    cf, time0, phases = calibrate_process2(loc, ids[0], initial)
    
    ## process all traces ##
    message("Begin processing 13C traces.", verbose)
    for i, item in enumerate(ids):
        message(f"Processing file {item}...", verbose)
        ppm, trace, timestamp = process_trace2(loc, item, cf, phases)
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
    writer = pd.ExcelWriter(loc + '/StackSpec.xlsx', engine='xlsxwriter')
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
    message(f"Completed successfully. Output: {loc}/StackSpec.xlsx", verbose)
    
def main():
    """Handle user input and call stack processing function."""
    print("Welcome to NMR processing pipeline.\n")
    correct = False
    while not correct:
        print("Please give the filepath of the experiment directory.")
        print("If the script was run from the experiment directory, " +
              "just press the [return] key")
        loc = input("")
        if loc == "":
            loc = os.getcwd()
        print("\nThe experiment directory is: " + loc)
        userinput = input("Is this correct? (y/n): ")
        if userinput == "y":
            correct = True
    correct = False
    while not correct:
        print("Please give the sequence of IDs for all " +
              "13C spectra to be processed.")
        print("For example: 12-52 even, 56, 61-67 odd, 70-73 all")
        userinput = input("")
        indices = []
        series = userinput.split(',')
        for serie in series:
            s = serie.strip()
            tokens = s.split(' ')
            if len(tokens) == 1:
                indices.append(tokens[0])
            else:
                bounds = tokens[0].split('-')
                if tokens[1] == 'all':
                    indices.extend(
                        [str(i) for i in range(int(bounds[0]), int(bounds[1]) + 1)])
                else:
                    indices.extend(
                        [str(i) for i in range(int(bounds[0]), int(bounds[1]) + 1, 2)])
        print("\nThe spectra to be processed are: " + ', '.join(indices))
        userinput = input("Is this correct? (y/n): ")
        if userinput == "y":
            correct = True
    init = input("\nPlease give the ID of the initial spectrum, " +
                 "to calibrate the timestamps: ")
    print("\nPlease review the processing parameters:")
    print("Directory: " + loc)
    print("13C spectrum IDs: " + ', '.join(indices))
    print("Initial timepoint reference: " + init)
    response = ""
    while response not in ['y', 'yes', 'no', 'n']:
        response = input("Are these parameters correct? (y/n): ")
    if response in ['y', 'yes']:
        process_stack(loc, indices, initial=init, verbose=True)
        return
    print("Please run the python script again and re-enter the parameters.")
    
if __name__ == "__main__":
    main()