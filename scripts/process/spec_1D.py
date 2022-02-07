#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 12:00:24 2021

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Process raw Bruker FID file into fourier transformed signal vector
   and display plot
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

import math, os, subprocess, sys, time
import nmrglue
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit as fit
from scipy.interpolate import UnivariateSpline

SCDIR = os.path.dirname(__file__)   # location of script

ver = "_2"

def process_spec(loc, item, isotope):
    """Get process parameters using nmrglue and write NMRpipe script

    Arguments:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor

    Returns: ppm correction to calibrate the x axis
    """
    runf = loc.split('/')[-1]

    calib = -1
    p0s = 999.0
    p1s = 999.0
    if isotope == '13C':
        ppm_max = 200.
        ppm_min = 0.
        if '13CGlc' in runf:
            refpk = 72.405 # glc
        elif '13CThr' in runf:
            refpk = 22.169
        elif '13CLeu' in runf:
            refpk = 23.596 # leu
        if runf.startswith('20210519'):
            p0s = 84.
            p1s = 81.
            calib = 69.950
        elif runf.startswith('20210525'):
            p0s = 81.
            p1s = 84.
            calib = 69.950        ## For 5/25 run
        elif runf.startswith('20210404'):
            p0s = 84.
            p1s = 81.
            calib = 69.355
        elif runf.startswith('20210322'):
            p0s = 84.
            p1s = 81.
            calib = 20.579
        elif runf.startswith('20210723'):
            p0s = 15.
            p1s = 100.
            calib = 69.884
        elif runf.startswith('20210731'):
            p0s = 84.
            p1s = 81.
            calib = 38.788
            refpk = 42.279
        elif runf.startswith('20210324'):
            p0s = 84.
            p1s = 81.
            calib = 19.079
        elif runf.startswith('20211203'):
            p0s = 84.
            p1s = 81.
            calib = 19.639
        # p0s = 81.0#76.0 #91.0 # 95
        # p1s = 84.0#84.0 # 68.0  #66

    elif isotope == '1H':
        if runf.startswith('20210324'):
            p0s = -135.
            p1s = 90.
            calib = 5.226
            refpk = 5.223 # acetate
        elif runf.startswith('20210519'):
            p0s = 128.
            p1s = 95.
            calib = 0.8594
            refpk = 0.8642
        elif runf.startswith('20211203'):
            p0s = 180.
            p1s = 10.
            calib = 5.216
            refpk = 5.223
        ppm_max = 12.
        ppm_min = 0.
    elif isotope == '31P':
        p0s = 15.0
        p1s = -117.0

    ## CONVERT and LOAD SPECTRUM ##
    subprocess.run(
        ["csh", SCDIR + "/bruker_" + isotope + ver + ".com", item, "f"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)

    ## LINE BROADENING, FT, BASELINE CORRECTION ##
    os.chdir(f"{loc}/{item}")
    time.sleep(1.)
    subprocess.run(["csh", SCDIR + "/ft_proc.com", item + "_f"])
    time.sleep(1.)
    os.chdir(loc)

    ## PHASE SHIFT ##
    dic, data = nmrglue.pipe.read(f"{loc}/{item}/{item}_f_1D.ft")

    #INSERT LOOP#
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
    #END INSERT#

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

    if calib < 0:
        uc = nmrglue.pipe.make_uc(dic, data, dim=0)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(uc.ppm_scale(), data.real, lw=0.5)
        ax.set_xlim(ppm_max, ppm_min)
        print("Find the ppm shift of the peak you wish to calibrate.")
        print("Note the x-coordinate, then close the window.")
        plt.show()
        calib = input("Type the ppm shift of the peak you wish to calibrate: ")
    # refpk = input("Type the reference ppm shift for the peak: ")
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
    subprocess.run(["csh", SCDIR + "/ft_proc.com", "13_f"])
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


def deconstruct(ppm, trace, pattern='15N'):
    """
    Deconstruct 13C splitting pattern of peak using various anchor points and
    preset splitting patterns.

    Parameters
    ----------
    ppm : list of floats
        The ppm shifts -- x data of spectrum.
    trace : list of floats
        The signal intensities -- y data of spectrum.
    anchor : float
        The ppm shift of the leftmost peak in the alanine splitting pattern.
    pattern : string
        The isotope splitting pattern to follow. Options: '13C', '15N'.

    Returns
    -------
    None, but displays overlay plot.

    """
    ## Display spec and find the location of leftmost peak
    ## This will be the anchor point of the remainder of the peak pattern
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(ppm, trace, lw=0.5)
    # ax.set_xlim(200.0, 0.0)
    # print("Find the ppm shift of the leftmost peak of Ala-C2.")
    # print("Note the x-coordinate, then close the window.")
    # plt.show()
    # anchor = input("Type the anchor ppm shift here: ")
    # pattern = input("Type the split type here (13C/15N): ")
    # anchor = float(anchor)
    anchor = 52.3425 # for 5/25 run

    idxk = [i for i in range(len(ppm)) if ppm[i]<anchor+0.2 and ppm[i]>anchor-0.8]
    ppm = [ppm[i] for i in idxk]
    trace = [trace[i] for i in idxk]
    ppm.reverse()
    trace.reverse()

    def logistic (L, k, x0, xarr):
        """Evaluate a logistic function given params."""
        try:
            result = np.asarray([L/(1. + math.exp(-1*k*(x - x0))) for x in xarr])
        except OverflowError:
            result = np.zeros_like(xarr)
        return result
    def gaussian(L, k, x0, xarr):
        """Evaluate a gaussian function given params."""
        try:
            result =  np.asarray([k*L*math.exp(-k*(x - x0))/(1. + math.exp(-1*k*(x - x0)))**2 for x in xarr])
        except OverflowError:
            result = np.zeros_like(xarr)
        return result
    def get_fcn(a, pattern, name='logistic'):
        c1 = -0.3588 # split distance from 13C-C1
        c3 = -0.2310 # split distance from 13C-C3
        h2 = 0.0 # split distance from C2-H
        n2 = -0.0371 # split distance from C2-N
        s23 = -0.1685 # shift from leftmost peak to left peak of 2-3 split
        c123 = [a, a + c1, a + c3, a + c1 + c3] # four-split
        c23  = [a, a + c3] # two-split
        c123_23 = [c123, [e + s23 for e in c23]]
        c123_23_n = [
            c123,
            [e + n2 for e in c123],
            [e + s23 for e in c23],
            [e + s23 + n2 for e in c23],
            [e + n2/2 for e in c123],
            [e + s23 + n2/2 for e in c23],
            ]
        c123_23_h = [
            c123,
            [e + h2 for e in c123],
            [e + s23 for e in c23],
            [e + s23 + h2 for e in c23]
            ]
        p = c123_23
        if pattern == '13C':
            p = c123_23
        elif pattern == '15N':
            p = c123_23_n
        print(p)
        fcn_superset = []
        if pattern == '13C':
            pass
        ftype = logistic
        if name == 'gaussian':
            ftype = gaussian
        if pattern == '15N':
            def fcn2(x, L1, L2, L3, L4, L5, L6, k1, k2, k3, k4, k5, k6):
                cpnts = []
                total = np.zeros_like(x)
                cpnts.append([ftype(L1, k1, h, x) for h in p[0]])
                cpnts.append([ftype(L2, k2, h, x) for h in p[1]])
                cpnts.append([ftype(L3, k3, h, x) for h in p[2]])
                cpnts.append([ftype(L4, k4, h, x) for h in p[3]])
                cpnts.append([ftype(L5, k5, h, x) for h in p[4]])
                cpnts.append([ftype(L6, k6, h, x) for h in p[5]])
                for ele in cpnts:
                    for arr in ele:
                        total = np.add(total, arr)
                return total
            def fcn3(L, k, i, x):
                total = np.zeros_like(x)
                result = [ftype(L, k, h, x) for h in p[i]]
                for arr in result:
                    total = np.add(total, arr)
                return total
            return fcn2, fcn3, p
        return None

    gfcn, gsfcns, gp = get_fcn(anchor, pattern, name='gaussian')
    cumtrace = np.insert(integrate.cumtrapz(trace, ppm), 0, 0.)
    gpopt, gpcov = fit(gfcn, ppm, trace, p0=[5000]*6+[150]*6)
    print(gpopt)
    gperr = np.sqrt(np.diag(gpcov))
    print(gperr)
    maxmaxcoeff = gpopt + gperr * np.asarray([1]*12)
    maxmincoeff = gpopt + gperr * np.asarray([1]*6 + [-1]*6)
    minmaxcoeff = gpopt + gperr * np.asarray([-1]*6 + [1]*6)
    minmincoeff = gpopt + gperr * np.asarray([-1]*12)
    maxarr = np.fmax(np.fmax(np.fmax(
        gfcn(ppm, *maxmaxcoeff),
        gfcn(ppm, *maxmincoeff)),
        gfcn(ppm, *minmaxcoeff)),
        gfcn(ppm, *minmincoeff))
    minarr = np.fmin(np.fmin(np.fmin(
        gfcn(ppm, *maxmaxcoeff),
        gfcn(ppm, *maxmincoeff)),
        gfcn(ppm, *minmaxcoeff)),
        gfcn(ppm, *minmincoeff))

    areas = []
    print("Calculating absolute areas:")
    for i in range(6):
        area = integrate.trapezoid(gsfcns(gpopt[i], gpopt[i+6], i, ppm), ppm)
        areas.append(area)
        print(f"  Peakset {i}: {area}")
    print("Calculating relative areas:")
    for i in range(6):
        print(f"  Peakset {i}: {areas[i]/sum(areas)}")
    print(f"(1,2,3-13C)-alanine: {(areas[0] + areas[1] + areas[4])/sum(areas)}")
    print(f"(2,3-13C)-alanine: {(areas[2] + areas[3] + areas[5])/sum(areas)}")
    print(f"(15N,13C)-alanine: {sum(areas[:4])/sum(areas)}")
    print(f"(14N,13C)-alanine: {sum(areas[4:])/sum(areas)}")
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    for l in gp:
        for i in l:
            ax1.axvline(i, ls=':')
    ax1.plot(ppm, trace, lw=2)
    ax1.plot(ppm, gfcn(ppm, *gpopt), lw=2)
    for i in range(6):
        ax1.plot(ppm, gsfcns(gpopt[i], gpopt[i+6], i, ppm), lw=2)
    ax1.plot(ppm, np.asarray(trace) - gfcn(ppm, *gpopt), lw=2)
    ax2 = fig.add_subplot(222)
    for l in gp:
        for i in l:
            ax2.axvline(i, ls=':')
    ax2.fill_between(ppm, minarr, maxarr, alpha=0.3, color="C2")
    ax2.plot(ppm, trace, lw=2)
    ax2.plot(ppm, gfcn(ppm, *gpopt), lw=2)
    ax3 = fig.add_subplot(223)
    for l in gp:
        for i in l:
            ax3.axvline(i, ls=':')
    ax3.plot(ppm, gfcn(ppm, *gpopt), lw=2, color='C1')
    for i in range(6):
        ax3.plot(ppm, gsfcns(gpopt[i], gpopt[i+6], i, ppm), lw=2, color=f"C{i+2}")
    ax4 = fig.add_subplot(224)
    for l in gp:
        for i in l:
            ax4.axvline(i, ls=':')
    ax4.plot(ppm, trace, lw=2)
    ax4.plot(ppm, gfcn(ppm, *gpopt), lw=2)
    ax4.plot(ppm, np.asarray(trace) - gfcn(ppm, *gpopt), lw=2, color="C8")

    ax1.set_title("Deconstructed and Raw NMR Signal")
    ax2.set_title("Curve-fit and Raw Signal w/ S.D.")
    ax3.set_title("Deconstructed Signal Components")
    ax4.set_title("Curve-fit, Raw, and Residual Signal")
    fig.supxlabel("¹³C chemical shift (ppm)")
    fig.supylabel("signal (unitless)")

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

    ppm, trace = process_spec(loc, fid, isotope)
    #deconstruct(ppm, trace)
    #subtract(loc, fid, isotope)

    #return

    #ppm, trace = process_spec(loc, fid, isotope)

    ## for Ala
    # trace = [trace[i] for i in range(len(trace)) if ppm[i] >= 52.425 and ppm[i] <= 54.425]
    # ppm = [ppm[i] for i in range(len(ppm)) if ppm[i] >= 52.425 and ppm[i] <= 54.425]
    ##
    lims = [200., 0.]
    if isotope == '1H':
        lims = [12., 0.]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ppm, trace, lw=0.5)
    ax.set_xlim(*lims) # 53.8, 53.0
    plt.show()

    """
    run_name = loc.split('/')[-1]
    with open(loc+'/'+run_name+'_alatrace.tsv', 'w') as wf:
        for i, p in enumerate(ppm):
            wf.write(f"{p}\t{trace[i]}\n")
    """

    print("Done.")
    return

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
