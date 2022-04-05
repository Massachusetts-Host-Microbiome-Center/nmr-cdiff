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

ver = '_2' # '_2'  # use '' or '_2'

def calibrate_process(loc, item, acq1, isotope, verbose=True):
    """Get process parameters using nmrglue and write NMRpipe script

    Arguments:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor

    Returns: ppm correction to calibrate the x axis
    """
    runf = loc.split('/')[-1]
    ## GET INITIAL TIMESTAMP ##
    time0 = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{acq1}/pulseprogram"))
    message("Loaded initial timestamp.", verbose)

    ## GET BRUKER PARAMS (NS) ##
    bruker_dic, _ = nmrglue.bruker.read(f"{loc}/{item}")
    n_scans = bruker_dic['acqus']['NS']

    ## CONVERT and LOAD SPECTRUM ##
    subprocess.run(
        ["csh", f"{SCDIR}/bruker_{isotope}{ver}.com", item, "calibr"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)
    message(f"Converted initial {isotope} spectrum for calibration.", verbose)
    calib = -1
    p0s = 999.0
    p1s = 999.0
    if isotope == '13C':
        ppm_max = 200.
        ppm_min = 0.
        if '13CGlc' in runf:
            refpk = 72.405
        elif '13CThr' in runf:
            refpk = 22.169
        if runf.startswith('20210519'):
            p0s = 84.
            p1s = 81.
            calib = 69.950
        elif runf.startswith('20210525'):
            p0s = 81.
            p1s = 84.
            calib = 69.950
        elif runf.startswith('20210404'):
            p0s = 84.
            p1s = 81.
            calib = 69.355
        elif runf.startswith('20210723'):
            p0s = 15.
            p1s = 100.
            calib = 69.884
        elif runf.startswith('20210324'):
            p0s = 84.
            p1s = 81.
            calib = 19.079
        elif runf.startswith('20210731'):
            p0s = 84.
            p1s = 81.
            calib = 38.788
            refpk = 42.279
        elif runf.startswith('20220325'):
            p0s = 84.
            p1s = 81.
            calib = 69.79
        # elif runf.startswith('20220321'):
        #     p0s = 84.
        #     p1s = 81.
            # calib = 69.79
        # p0s = 81.0#76.0 #91.0 # 95
        # p1s = 84.0#84.0 # 68.0  #66
    elif isotope == '1H':
        if runf.startswith('20210324'):
            p0s = -135.
            p1s = 90.
            calib = 5.226
        refpk = 5.223
        ppm_max = 12.
        ppm_min = 0.
    elif isotope == '31P':
        p0s = 15.0
        p1s = -117.0

    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/ft_proc.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft")
    message("Completed LB, FT, BL transformations.", verbose)

    ## CALIBRATE PHASE CORRECTION IF NECESSARY ##
    if p0s == 999.0 and p1s == 999.0:
        p0t = p0s
        p1t = p1s
        while True:
            print(p0t, p1t)
            dic1, data1 = nmrglue.process.pipe_proc.ps(dic, data, p0=p0t, p1=p1t)
            uc = nmrglue.pipe.make_uc(dic1, data1, dim=0)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(uc.ppm_scale(), data1.real, lw=0.5)
            ax.set_xlim(ppm_max, ppm_min)
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
        p0s = p0t
        p1s = p1t

    ## PHASE SHIFT ##
    dic, data = nmrglue.process.pipe_proc.ps(dic, data, p0=p0s, p1=p1s)
    nmrglue.pipe.write(f"{loc}/{item}/{item}_calibr_1D.ft.ps", dic, data, overwrite=True)
    message("Completed phase correction.", verbose)

    ## BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/bl_proc.com", item + "_calibr"])
    time.sleep(1.)
    os.chdir(loc)
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_calibr_1D.ft.ps.bl")
    message("Completed baseline correction.", verbose)

    ## CALIBRATE PPM SHIFT IF NECESSARY ##
    if calib < 0:
        uc = nmrglue.pipe.make_uc(dic, data, dim=0)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(uc.ppm_scale(), data.real*(n_scans**0.5), lw=0.5)
        ax.set_xlim(ppm_max, ppm_min)
        print("Find the actual ppm shift of the glucose peak with reference shift 72.405.")
        print("Note the x-coordinate, then close the window.")
        plt.show()
        calib = input("Type the ppm shift of the peak you wish to calibrate: ")

    cf = float(refpk) - float(calib)
    message("Corrected PPM shift.", verbose)

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

    ## GET BRUKER PARAMS (NS) ##
    bruker_dic, _ = nmrglue.bruker.read(f"{loc}/{item}")
    n_scans = bruker_dic['acqus']['NS']

    ## CONVERT BRUK to PIPE ##
    subprocess.run(["csh", f"{SCDIR}/bruker_{isotope}{ver}.com", item],
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

    ## SHIFT AND NORMALIZE DATA ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_1D.ft.ps.bl")
    uc = nmrglue.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale() + cf
    trace = data.real * (n_scans**0.5)

    return ppm, trace, tim

def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, isotope='13C', c_each=False, verbose=True):
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
    cf, time0, phases = calibrate_process(loc, ids[0], initial, isotope, verbose)

    ## process all traces ##
    message("Begin processing " + isotope + " traces.", verbose)
    for i, item in enumerate(ids):
        message(f"Processing file {item}...", verbose)
        if c_each:
            cf, _, phases = calibrate_process(loc, item, initial, isotope, verbose)
        ppm, trace, timestamp = process_trace(loc, item, cf, phases, isotope)
            # if i == 0:
            #     stack.append(ppm)
            # stack.append(trace)
        if i == 0:
            stack.append([shift for shift in ppm if shift <= 200 and shift >= 0])
        start_index = next(i for i, shift in enumerate(ppm) if shift < 200)
        idx = (start_index, start_index+len(stack[0]))
        stack.append(trace[idx[0]:idx[1]])
        delta_t = (timestamp - time0).total_seconds()/3600.
        header[0].append(item)
        header[1].append(delta_t)
        it = iter(stack)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            print("PPM mismatch. Abort.")
            return

    stem = loc.split('/')[-1]

    ## write excel file ##
    message("Trace processing complete.", verbose)
    message("Writing full stack to Excel file.", verbose)
    stack_array = np.array(stack)
    df = pd.DataFrame(stack_array[1:], columns=stack_array[0], index=header[1])
    writer = pd.ExcelWriter(f"{loc}/{stem}_{isotope}.xlsx", engine='xlsxwriter')
    df = df.T
    # df['mean'] = df.mean(axis=1)
    # header[0].append('mean')
    # header[1].append('mean')
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
    message(f"Completed successfully. Output: {loc}/{stem}_{isotope}.xlsx", verbose)

def detect_spectra(isotope='13C', init=11):
    """Returns a list of spectra IDs for the given isotope."""
    iso_spectra = []
    dirs = glob.glob('*/')
    st = [str(y) for y in sorted([int(x[:-1]) for x in dirs])]
    for fn in st:
        if int(fn) >= int(init):
            try:
                dic, data = nmrglue.bruker.read(fn)
                udic = nmrglue.bruker.guess_udic(dic, data)
                if udic[0]['label'] == isotope and udic[0]['encoding'] == 'direct':
                    iso_spectra.append(fn)
            except OSError:
                pass
    return iso_spectra

def main(iso, init, calibrate_each):
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
    process_stack(loc, indices, initial=init, isotope=iso,
            c_each=calibrate_each, verbose=True)

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
        if args[2] == '-ce':
            return iso, init, True
        print("Ignoring excess arguments: " + " ".join(args[2:]))
    return iso, init, False

if __name__ == "__main__":
    main(*handle_args(sys.argv[1:]))
