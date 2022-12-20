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

import collections
import copy
import datetime
import glob
import os
import subprocess
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np
import pandas as pd
import scipy
# import skimage

SCDIR = os.path.dirname(__file__)   # location of script
CURVE_FIT_RESOLUTION = 2**19

# ver = ''
# ver = '_2'
ver = '_3'
# ver = '_4'

ISOTOPE_PARAMS = {
    '1H': {
        'fit_ppm_delta': 0.01,   # +/- drift tolerance for peak centers during fitting
        'fit_fwhm_max': 0.25,   # max FWHM of fit curve
        'fit_fwhm_min': 0.001,  # min FWHM of fit curve (unused, no lower bound)
        'plot_bounds': (12., 0.),   # PPM shift bounds for plotting
        'assignment_display_window': 1.,    # extra plot padding to contextualize peaks during assignment
        'fit_cluster_msep': 0.1,    # minimum separation for distinct clusters, also cluster buffer for automatic peak assignment
        'fit_cluster_reach': 0.05,  # buffer to include around cluster for curve fitting
    },
    '13C': {
        'fit_ppm_delta': 5,
        'fit_fwhm_max': 0.5,
        'fit_fwhm_min': 0.01,
        'plot_bounds': (200., 0.),
        'assignment_display_window': 8,
        'fit_cluster_msep': 10,
        'fit_cluster_reach': 8,
    }
}

def plot_datasets(xx, yys, ppm_bounds=(200., 0.)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for yy in yys:
        ax.plot(xx, yy, lw=0.5)
    ax.set_xlim(*ppm_bounds)
    plt.show()

def plot_with_labels(xx, yy, labeli, labels, ppm_bounds=(200., 0.), yopt=None, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xx, yy, lw=0.5)
    if yopt is None:
        yopt = yy[labeli]
    ax.plot(xx[labeli], yopt, ls='', marker='*', **kwargs)
    for ci, xi, yi in zip(labels, xx[labeli], yopt):
        ax.annotate(ci, (xi,yi), textcoords='offset points', xytext=(0,10),
            ha='center')
    ax.set_xlim(*ppm_bounds)
    plt.show()

def to_stream(dic, data):
    """Prepare data for NMRpipe.
    Adapted from nmrglue.fileio.pipe
    """
    if data.dtype == "complex64":
        data = ng.fileio.pipe.append_data(data)
    data = ng.fileio.pipe.unshape_data(data)

    # create the fdata array
    fdata = ng.fileio.pipe.dic2fdata(dic)

    # write to pipe
    if data.dtype != 'float32':
        raise TypeError('data.dtype is not float32')
    if fdata.dtype != 'float32':
        raise TypeError('fdata.dtype is not float32')

    datastream = fdata.tobytes() + data.tobytes()
    return datastream

def get_ppm_bounds(isotope):
    if isotope == '13C':
        return (200., 0)
    else:
        return (12., 0.)

def get_timestamp(fidpath):
    """Get FID timestamp from acqus file."""
    dic, _ = ng.bruker.read(fidpath)
    return datetime.datetime.fromtimestamp(dic['acqus']['DATE'])

def pipe_process(isotope, fid):
    """Process Bruker FID using NMRPipe and load using nmrglue."""
    pipe_output = subprocess.run(
        ["csh", f"{SCDIR}/proc_{isotope}{ver}.com", fid],
        stdout=subprocess.PIPE)
    return ng.fileio.pipe.read(pipe_output.stdout)

def pipe_bl(dic, data):
    """Perform baseline correction using NMRPipe."""
    pipe_output = subprocess.run(
        ["csh", f"{SCDIR}/bl.com"],
        input=to_stream(dic, data),
        stdout=subprocess.PIPE)
    return ng.fileio.pipe.read(pipe_output.stdout)

def calibrate_process(loc, item, acq1, isotope, verbose=True):
    """Get process parameters using NMRglue and write NMRpipe script

    Parameters:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor

    Returns: ppm correction to calibrate the x axis
    """
    
    ppm_bounds = get_ppm_bounds(isotope)

    runf = loc.split('/')[-1]

    dic, _ = ng.bruker.read(f'{loc}/{acq1}') # get initial timestamp
    time0 = datetime.datetime.fromtimestamp(dic['acqus']['DATE'])

    dic, data = ng.bruker.read(f'{loc}/{item}') # get number of scans
    n_scans = dic['acqus']['NS']

    dic, data = pipe_process(isotope, item) # CONVERT, PROCESS, and LOAD SPECTRUM

    p0, p1 = ng.process.proc_autophase.manual_ps(data, notebook=False) # MANUAL PHASE SHIFT
    dic, data = ng.process.pipe_proc.ps(dic, data, p0=p0, p1=p1)

    dic, data = ng.process.pipe_proc.mult(dic, data, r=n_scans, inv=True) # NORMALIZE BY NUMBER OF SCANS
    
    dic, data = pipe_bl(dic, data) # BASELINE CORRECTION
    # if dic['FDF2LABEL'] == "1H":
    #     data = ng.process.proc_bl.baseline_corrector(data, wd=25)

    ## CALIBRATE PPM SHIFT ##
    print("Find the experimental ppm shift of a reference peak to " \
            "calibrate the chemical shift axis.")
    ref_shift = input("Type the reference ppm shift of the peak to " \
        "calibrate: ")
    ref_shift = float(ref_shift)
    print("Now, find the experimental shift of that peak.")
    print("Note the x-coordinate, then close the window.")
    uc = ng.pipe.make_uc(dic, data, dim=0)
    plot_datasets(uc.ppm_scale(), [data.real], ppm_bounds=ppm_bounds)
    calib = input("Type the ppm shift of the peak you wish to calibrate: ")
    actual_shift = float(calib)

    calibration_shift = float(ref_shift) - float(actual_shift)
    print(f"Shifting by {calibration_shift} ppm.")
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    return calibration_shift, time0, [p0, p1], ppm

def rmse(vector):
    """Calculate the Root Mean Squared Error (RMSE) of a 1D-array."""
    return np.sqrt(((vector - vector.mean()) ** 2).mean())

def ppm_ivl(uc, xx): 
    """Convert PPM scale to index."""
    return abs(uc.f(0, unit='ppm') - uc.f(xx, unit='ppm'))

def shift_params(s, *args):
    """Shift parameters by a factor s.
    Only acts on first element, which is usually the ppm shift.
    """
    for params in args:
        for param in params:
            param[0] -= s
    return args

def sg_peakpick(data, w_len, p_order, rsg=2, rd=2, nr=None, plot=True):
    """Detect peaks using a second-derivative Savitzky-Golay filter.
    
    TODO: calculate SD from noise region, not whole spec

    Parameters:
    data -- FT-NMR spectrum for peak-picking.
    w_len -- SG filter window size, must be odd. Larger window size results in greater
             smoothing.
    p_order -- order of SG polynomial.
    rsg -- number of standard deviations for peak-finding minimum prominence, SG.
    rd -- number of standard deviations for peak-finding minimum prominence, FT data.
    nr -- two-ple of index bounds of noise region for std calculation.
    
    Return lost containing indices of found peaks and estimated widths.
    """
    # Compute Savitsky-Golay filter, 2nd derivative
    fd2 = -1*scipy.signal.savgol_filter(data, w_len, p_order, deriv=2)

    # Calculate noise thresholds
    if nr is None:
        err_f = rd*np.std(data)
        err_fd2 = rsg*np.std(fd2)
    else:
        err_f = rd*np.std(data[slice(*nr)])
        err_fd2 = rsg*np.std(fd2[slice(*nr)])

    # Find peaks, estimate widths and amplitudes, filter by noise thresholds
    pks, prp = scipy.signal.find_peaks(fd2, distance=3, width=(3, 75), prominence=err_fd2)
    pkw = np.array([abs(prp['right_ips'][i] - prp['left_ips'][i]) for i in range(pks.shape[0])])
    amps = np.array(data[[int(pk) for pk in pks]])
    pks = pks[amps > err_f]
    pkw = pkw[amps > err_f]
    x = np.array(range(len(data), 0, -1))

    # Plot found peaks
    if plot:
        pklabels = list(range(pks.shape[0]))
        plot_with_labels(x, fd2, pks, pklabels, ppm_bounds=(max(x), min(x)))

    return list(pks), list(pkw)

def peak_pick(dic, data, plot=True):
    """Use Savitzky-Golay method to detect both fine and coarse peaks.
    
    TODO: tune for 13C
    
    Parameters:
    data -- vector of data for peak-picking.
    isotope -- isotope to determine SG filter parameters.

    Returns list containing positions of peaks detected.
    """
    isotope = dic['FDF2LABEL']
    uc = ng.pipe.make_uc(dic, data, dim=0)

    # Perform peak-peaking using Savitsky-Golay filter tuned to coarse and fine detail
    if isotope == '1H':
        nr = (int(uc.f(12, unit="ppm")), int(uc.f(8, unit="ppm")))
        coarse, cpeakw = sg_peakpick(data, 41, 2, rsg=5, rd=4, nr=nr, plot=plot)
        fine, fpeakw = sg_peakpick(data, 11, 2, rsg=8, rd=8, nr=nr, plot=plot)
    # elif isotope == '13C': # TODO: OPTIMIZE FOR 13C
    #     nr = (int(uc.f(160, unit="ppm")), int(uc.f(130, unit="ppm")))
    #     coarse_peaks = list(sg_peakpick(data, 31, 3, r=2, nr=nr))
    #     fine_peaks = list(sg_peakpick(data, 11, 2, r=5, nr=nr))
    all_peaks = fine + coarse   # ppm shifts of peaks
    all_widths = fpeakw + cpeakw    # estimated widths of peaks
    peaks = []
    peakw = []
    for pair in zip(all_peaks, all_widths):
        peaks.append((pair[0],))
        peakw.append((pair[1],))
    return peaks, peakw

def peak_fit(dic, data, r=6, sep=0.005, plot=True, vb=False):
    """Peak-pick NMR spectrum.
    Includes peak-picking, peak assignment, and integration functionality. The
    reference peaks to be fit should be specified in cfg_{isotope}.txt and
    provided in the directory of the run to be processed.

    Parameters:
    ppm -- chemical shift axis of the data, in ppm
    data -- signal axis of the data
    r -- RMSE noise threshold for peak detection (default: 6)
    plot -- True (default) to plot fit peaks with cluster labels over the data.
    vb -- True for verbose mode.

    Returns a dictionary mapping compounds to their peak areas, and a curve
    summing all of the fit peaks.
    """
    # ----------------------------------- UTILITIES ------------------------------------
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    data = data.real
    isotope = dic["FDF2LABEL"]
    plot_bounds = ISOTOPE_PARAMS[isotope]['plot_bounds']
    if isotope == '1H':
        data = -1*data

    def p2f_interval(xx): return abs(uc.f(0, unit='ppm') - uc.f(xx, unit='ppm')) # convert PPM scale to index
    def p2i_interval(xx): return int(p2f_interval(xx)) # convert PPM scale to index # convert PPM scale to index
    def p2i(xx): return uc.f(xx, unit='ppm')
    def i2p(yy): return uc.ppm(yy) # convert index to PPM shift

    # ---------------------------------- PEAK-PICKING ----------------------------------
    if isotope == '1H':
        pthres = 1000 # 300
        shifts, peakw = peak_pick(dic, data, plot=plot)
        cluster_sep = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_cluster_msep'])
        distances = scipy.spatial.distance.cdist(shifts, shifts)
        nclust, cIDs = scipy.sparse.csgraph.connected_components(
            (distances <= cluster_sep).astype(int), 
            directed=False
        )
    else:
        pthres = r*rmse(data[(ppm <= 160) & (ppm >= 130)])
        shifts, cIDs, params, amps = ng.analysis.peakpick.pick( # Find peaks
            data, pthres=pthres, msep=(p2i_interval(sep),),
            algorithm='thres-fast', est_params=True, lineshapes=['g'],
            cluster=True, c_ndil=p2i_interval(0.3), table=False
        )

    if len(shifts) == 0:
        # Exit if no peaks detected
        print("No peaks detected; exit.")
        return {}, []    
    amps = np.array([data.real[a] for a in shifts])
    shifts_tup = np.array(shifts, dtype=np.int32)
    shifts = np.array([sh[0] for sh in shifts])
    peakw = np.array(peakw)

    if plot:
        plot_with_labels(ppm, data, shifts, cIDs, ppm_bounds=plot_bounds)
        plt.show()

    # --------------------------- ASSIGN PEAKS TO COMPOUNDS ----------------------------
    # Prepare arrays and constants for use within the loop
    refsh = pd.read_csv(f"cfg_{isotope}.txt", sep='\t', names=['Shift', 'Compound']) # Shift | Compound
    refpks = np.array([p2i(v) for _, v in refsh['Shift'].items()])
    assignments = collections.defaultdict(list) # dict to store compound assignments
    display_window = p2f_interval(ISOTOPE_PARAMS[isotope]['assignment_display_window'])

    # Assign compound to peaks by proximity to reference shifts, handling collisions
    # Iterate: cluster ID
    for cid in np.unique(cIDs):
        subpeaks = shifts[cIDs == cid] # detected peaks w/in cluster
        clb = subpeaks.min() - cluster_sep # calculate cluster bounds
        cub = subpeaks.max() + cluster_sep

        # Determine which reference peaks fall within cluster bounds
        close_refpks = [i for i, xx in enumerate(refpks) if xx > clb and xx < cub]

        # Assign peaks to compounds; skip clusters with no nearby reference shifts
        if len(close_refpks) == 1:
            # If only one matched ref. peak, assign all cluster sub-peaks to that compound
            assignments[refsh.at[close_refpks[0], 'Compound']].extend(subpeaks)
        elif len(close_refpks) > 1: # Collision!
            # If multiple ref. peaks, manually assign sub-peaks in plot window
            colliding_refs = refsh.loc[close_refpks, :]
            print(f"Cluster {cid} assigned to {len(close_refpks)} conflicting reference peaks.")
            print("In the plot, find the contributions from each " \
                    "peak using the reference shifts listed:")
            with pd.option_context(
                    'display.max_rows', None, 'display.max_columns', None):
                print(colliding_refs)
            for cpd in np.unique(colliding_refs['Compound']):
                while True:
                    plot_with_labels(
                        ppm, data.real, subpeaks, range(len(subpeaks)),
                        ppm_bounds=(i2p(min(subpeaks) - display_window),
                                    i2p(max(subpeaks) + display_window))
                    )
                    raw_text = input(f"Subpeak contributions of {cpd}: ")
                    try:
                        indices = [int(i) for i in raw_text.strip().split()]
                        break
                    except ValueError:
                        print("Invalid assignments.")
                if len(indices) > 0:
                    assignments[cpd].extend(subpeaks[indices])

    if len(assignments) == 0:
        print("No reference peaks detected; exit.")
        return {}, []

    # Filter out excluded data and prepare data structures for curve-fitting
    all_cpds = []
    used_pks = []
    for cpd, apks in assignments.items():
        if len(apks) > 0:
            all_cpds.append(cpd)
            used_pks.extend(apks)
    all_cpds = np.array(all_cpds)
    mask = np.isin(shifts, used_pks)
    shifts = shifts[mask]
    peakw = np.array(peakw)[mask]
    shifts_tup = shifts_tup[mask]
    cIDs = cIDs[mask]

    # --------------------------- CURVE-FIT PEAKS IN WINDOWS ---------------------------
    # Retreive reasonable curve-fitting parameter bounds for isotope
    ppm_delta = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_ppm_delta'])
    fwhm_min = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_fwhm_min'])
    fwhm_max = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_fwhm_max'])
    amps = np.array([data.real[int(a)] for a in shifts], dtype=np.float64)
    amp_bounds = np.array([(0.05*b, 1.2*b) for b in amps])

    # Define start parameters and bounds for each peak
    cv = 'l'
    if cv == 'v':
        cparams = np.array([((a, b[0], b[0]),) for a, b in zip(shifts, peakw)])
        cbounds = []
        for par in cparams:
            cbounds.append(((
                (par[0][0] - ppm_delta, par[0][0] + ppm_delta),    # shift bounds
                (None, fwhm_max), # gauss bounds
                (None, fwhm_max), # lorentz bounds
            ),))
    elif cv == 'l':
        cparams = np.array([((a, 2*b[0]),) for a, b in zip(shifts, peakw)])
        cbounds = []
        for par in cparams:
            cbounds.append(((
                (par[0][0] - ppm_delta, par[0][0] + ppm_delta),    # shift bounds
                (None, p2f_interval(fwhm_max)), # lorentz bounds
            ),))
    cbounds = np.array(cbounds)

    # Prepare arrays and constants for use within the loop
    curve_params = [cparams, amps, cbounds, amp_bounds, shifts_tup, cIDs]
    fit_params = np.zeros_like(cparams) # array for fit result, curve shift + FWHM
    fit_amps = np.zeros_like(amps)      # array for fit result, curve amplitudes
    cluster_sep = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_cluster_msep'])
    cluster_reach = p2f_interval(ISOTOPE_PARAMS[isotope]['fit_cluster_reach'])

    # Calculate curve fitting windows and iterate for curve fit
    distances = scipy.spatial.distance.cdist(shifts_tup, shifts_tup)
    n_fitwd, fitwd = scipy.sparse.csgraph.connected_components(
        (distances <= cluster_sep).astype(int), 
        directed=False
    )
    for wid in range(n_fitwd):
        mask = (fitwd == wid)
        fit_args = [pset[mask] for pset in curve_params]
        
        # slice region from +/- <cluster_reach> for curve-fitting
        lb = min(shifts[mask]) - cluster_reach
        ub = max(shifts[mask]) + cluster_reach
        _, ndata = ng.process.pipe_proc.ext(dic, data, x1=lb, xn=ub)

        # Shift peak parameters into slice and perform curve fit
        fit_args[4] = shift_params(lb, fit_args[4])[0]  # Shift ppm shift
        fit_args[0] = shift_params(lb, *fit_args[0])    # Shift VP params
        fit_args[2] = shift_params(lb, *fit_args[2])    # Shift VP param bounds
        nparams, namps, _ = ng.analysis.linesh.fit_spectrum(
            ndata, 
            [cv],
            *fit_args,
            p2i_interval(ppm_delta), 
            False, 
            verb=False,
            maxfev=1000000,
        )
        nparams = shift_params(-1*lb, *np.array(nparams)) # rescale result to full PPM axis

        # Deposit curve-fit parameters
        fit_params[mask] = np.array(nparams)
        fit_amps[mask] = np.array(namps)

    ## Prepare constants and arrays for plotting and integration
    resolution = CURVE_FIT_RESOLUTION
    scale = resolution/data.shape[0]
    sppm = np.linspace(ppm[0], ppm[-1], num=resolution)
    rppm = np.linspace(ppm[0], ppm[-1], num=data.shape[0])
    residuals = np.copy(data)

    # -------------- INTEGRATE, PLOT FIT PEAKS WITH COMPOUND ASSIGNMENTS ---------------
    signals = {}
    curve = np.zeros((resolution,))
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5, zorder=2)

    # Simulate, integrate, and plot curve for each compound
    for cpd in all_cpds:
        # Select peaks belonging to compound
        include = np.isin(shifts, assignments[cpd])
        idxs = shifts[include]

        # Simulate spectrum with fit parameters
        sim_smooth = ng.analysis.linesh.sim_NDregion( # Simulate spectrum in high-resolution
            (resolution,), 
            [cv], 
            scale*fit_params[include], 
            fit_amps[include],
        )

        # Integrate curve, store result, and add to overall signal
        areacurve = -1*np.trapz(sim_smooth, x=sppm)
        signals[cpd] = areacurve
        curve += sim_smooth

        # Plot simulated curve and annotated peak summits; update residual
        if plot:
            color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(sppm, sim_smooth, lw=1, color=color, zorder=4)
            ax.plot(ppm[idxs], fit_amps[include], ls='', marker='*', color=color, zorder=5)
            for xx, yy in zip(ppm[idxs], fit_amps[include]):
                ax.annotate(cpd, (xx,yy), textcoords='offset points',
                    xytext=(0,10), ha='center', zorder=6)
            sim_native = ng.analysis.linesh.sim_NDregion( # Simulate spectrum in native resolution
                (data.shape[0],), 
                [cv], 
                fit_params[include], 
                fit_amps[include],
            )
            residuals -= sim_native

    # Plot simulated signal + residual and show plot
    if plot:
        ax.plot(sppm, curve, lw=0.5, color='k', zorder=3)
        ax.plot(rppm, residuals, lw=0.5, color='gray', zorder=1)
        ax.set_xlim(*plot_bounds)
        plt.show()

    # Return curve integrations and simulated spectrum
    return signals, curve

def ridge_trace(dic, data, wr=0.01, plot=True):
    """Simple function to get the heights of reference peaks in a spectrum.
    Takes the maximum signal within window of each peak. Does not functionally
    trace ridge across spectra; cannot handle drifting chemical shifts, eg. due
    to pH changes. Shift must remain within window of reference throughout time
    course.

    Parameters:
    dic -- "universal" dictionary for spectrum
    data -- data vector for spectrum
    wr -- window radius, distance from reference shift to include in window

    Returns a dictionary containing the amplitude of each peak.
    """
    uc = ng.pipe.make_uc(dic, data, dim=0)
    isotope = dic["FDF2LABEL"]
    def p2i(xx): return int(uc.f(xx, unit='ppm'))

    signals = {}
    shifts = []
    amps = []
    cpds = []
    with open(f"cfg_{isotope}.txt", 'r') as rf:
        for line in rf:
            fields = line.strip().split('\t', 1)
            shift = float(fields[0])
            cpd = fields[1]
            bounds = (p2i(shift-wr), p2i(shift+wr))
            peak = f"{fields[1]} {shift}"
            amplitude = max(data[bounds[1]:bounds[0]])
            signals[peak] = amplitude
            shifts.append(shift)
            cpds.append(cpd)
            amps.append(amplitude)
    if plot:
        plot_bounds = ISOTOPE_PARAMS[isotope]['plot_bounds']
        plot_with_labels(uc.ppm_scale(), data, [p2i(s) for s in shifts],
            [str(s) for s in shifts], ppm_bounds=plot_bounds, yopt=amps)
        plot_with_labels(uc.ppm_scale(), data, [p2i(s) for s in shifts],
            [c for c in cpds], ppm_bounds=plot_bounds, yopt=amps, lw=2)

    return signals

def process_trace(loc, item, cf, ps, overwrite=True, man_ps=False):
    """Process a single 1D 13C-NMR trace.

    Parameters:
    loc -- path to folder containing traces
    item -- ID of trace to process
    cf -- correction factor for ppm shift
    ps -- list with the default p0, p1 values for phase correction
    overwrite -- whether to process each spectrum. If False, reads the existing
        "ft" file that has already been processed.
    man_ps -- whether to prompt for manual phase correction

    Returns:
    ppm --  list of ppm values
    trace -- list of signal intensities
    timestamp --  timestamp of trace, calculated from acquisition data
    signals -- integrated signals (or raw signals) of compounds calculated from peak fit
    n_scans -- number of scans for spectrum
    curve -- simulated spectrum (or raw spectrum)
    """

    dic, _ = ng.bruker.read(f"{loc}/{item}")  # GET BRUKER PARAMS
    n_scans = dic['acqus']['NS']
    isotope = dic['acqus']['NUC1']
    timestamp = datetime.datetime.fromtimestamp(dic['acqus']['DATE'])

    if os.path.exists(f"{loc}/{item}/ft") and not overwrite:
        dic, data = ng.fileio.pipe.read(f"{loc}/{item}/ft") # Load a previously processed spectrum
    else:
        dic, data = pipe_process(isotope, item)  # Convert, LB, ZF, FT

        dic, data = ng.process.pipe_proc.ps(dic, data, p0=ps[0], p1=ps[1]) # Coarse PS and manual adjustment
        if sum(data) < 0:
            dic, data = ng.process.pipe_proc.ps(dic, data, p0=180, p1=0)
        if man_ps:
            p0, p1 = ng.process.proc_autophase.manual_ps(data, notebook=False)
            dic, data = ng.process.pipe_proc.ps(dic, data, p0=p0, p1=p1)

        dic, data = ng.process.pipe_proc.mult(dic, data, r=n_scans, inv=True) # Normalize by number of scans

        # dic, data = pipe_bl(dic, data)  # Baseline correction
        # if isotope == "1H":
        #     data = ng.process.proc_bl.baseline_corrector(data, wd=25)
        #     data = np.float32(data)

        # Calibrate PPM shift
        ppm = ng.pipe.make_uc(dic, data, dim=0).ppm_scale()
        if cf != 0:
            shift = cf*data.shape[0]/(max(ppm) - min(ppm))
            dic, data = ng.process.pipe_proc.ls(dic, data, shift, sw=False)

        # Write for future use
        ng.fileio.pipe.write(f"{loc}/{item}/ft", dic, data, overwrite=True)

    # Pick, fit, and integrate peaks
    ppm = ng.pipe.make_uc(dic, data, dim=0).ppm_scale()
    if isotope == "13C":
        signals, curve = peak_fit(dic, data, r=4, plot=False, vb=True)
        # signals = ridge_trace(dic, data, plot=True)
        # curve = data.real
    else:
        # Flatten water resonances
        # ppm_array = np.array(ppm)
        data[(ppm > 4.2) & (ppm < 4.9)] = 0

        # Pick peaks and calculate signal
        signals, curve = peak_fit(dic, data, r=12, sep=30/2*0.005, plot=True, vb=True) # r=4 # sep=30/2*0.001
        # signals = ridge_trace(dic, data, plot=True)
        curve = data.real

    return ppm, data.real, timestamp, signals, n_scans, curve

def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, iso='13C', vb=True, dry_run=False):
    """Process entire NMR stack and write to xlsx file.

    Parameters:
    loc -- path to folder containing traces
    ids -- list of IDs
    initial -- ID of initial acquisition, if different from first ID \
        (to find initial timepoint)
    vb -- if True, print status messages
    dry_run -- if True, does not write the Excel file
    """
    message("Begin NMR stack process.", vb)
    if initial == None:
        initial = ids[0]

    header = [[], []]
    stack = []
    curves = []

    ## calibrate and get correction factor ##
    message("Calibrate autophase and PPM shift...", vb)
    cf, time0, phases, c_ppm = calibrate_process(loc, ids[0], initial, iso, vb)
    c_ppm = c_ppm[(c_ppm <= 200) & (c_ppm >= 0)]
    stack.append(c_ppm)

    ## process all traces ##
    message("Begin processing " + iso + " traces.", vb)
    for i, item in enumerate(ids):
        message(f"Processing file {item}...", vb)
        # Process spec and chop signal vector
        ppm, trace, ts, signals, n_scans, curve = process_trace(loc, item, cf, phases)
        start_index = next(i for i, shift in enumerate(ppm) if shift < 200)
        idx = (start_index, start_index+len(stack[0]))
        vector = trace[idx[0]:idx[1]]
        curves.append(curve)

        # Store data and check for continuity
        stack.append(vector)
        delta_t = (ts - time0).total_seconds()/3600.
        header[0].append(item)
        header[1].append(delta_t)
        if i == 0:
            traj = pd.DataFrame(
                np.zeros((len(ids), len(signals)+2)),
                columns=['Time', 'Scans'] + [cpd for cpd in signals])
        signals['Time'] = delta_t
        signals['Scans'] = n_scans
        for cpd, cnc in signals.items():
            traj.at[i, cpd] = cnc
        it = iter(stack)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            print("PPM mismatch. Abort.")
            return

    stem = loc.split('/')[-1]

    fig = plt.figure(figsize=(15,11))
    ax = plt.axes(projection='3d')
    ax.grid(False)
    if iso == "13C":
        xrange = ppm.shape[0]
        xmin = 0
        xmax = 100
        lcol = 'indigo'
    else:
        xrange = ppm.shape[0]
        xmin = 0
        # xmax = 6
        xmax = 4.4
        lcol = 'turquoise'
    pp_array = np.linspace(ppm[0], ppm[-1], num=xrange)
    plot_ppm = pp_array[(pp_array <= xmax) & (pp_array >= xmin)]
    verts = []
    # x = plot_ppm
    # y = header[1]
    # X, Y = np.meshgrid(x, y)
    # print(X)
    # print(Y)
    # Z = np.array([cv[(pp_array <= xmax) & (pp_array >= xmin)] for cv in curves])
    # zmax = max(Z.max(), -1*Z.min())
    # norm = mpl.colors.CenteredNorm(halfrange=1000)
    # surf = ax.plot_surface(X, Y, Z, cmap=mpl.cm.coolwarm_r, norm=norm, antialiased=True)
    for i in range(len(curves)-1, -1, -1):
        curve = curves[i][(pp_array <= xmax) & (pp_array >= xmin)]
        ax.plot3D(plot_ppm, [header[1][i] for t in plot_ppm], curve, lw=0.3)
        verts.append(list(zip(plot_ppm, curve)))
    poly = mpl.collections.PolyCollection(verts, edgecolors=lcol, linewidths=0.3)
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=header[1][::-1], zdir='y')
    ax.set_xlabel(f"{iso} chemical shift (ppm)")
    ax.set_xlim3d(xmax, xmin)
    ax.set_xlim3d(xmin, xmax)
    ax.set_ylabel("Time (h)")
    ax.set_ylim3d(0, max(header[1]))
    ax.set_zlabel("NMR signal (unitless)")
    ax.set_zlim3d(min([min(cv) for cv in curves]), max([max(cv) for cv in curves]))
    ax.view_init(elev=-10., azim=-82.)
    ax.set_proj_type('ortho')
    ax.set_axis_off()
    plt.show()
    plt.close('all')

    ax = traj.plot(x='Time', y=list(traj.columns)[2:], marker='.')
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("NMR signal (unitless)")
    plt.show()
    plt.close('all')

    ## write excel file ##
    message("Trace processing complete.", vb)
    if not dry_run:
        message("Writing full stack to Excel file.", vb)
        stack_array = np.array(stack)
        df = pd.DataFrame(stack_array[1:], columns=stack_array[0], index=header[1])
        writer = pd.ExcelWriter(f"{loc}/{stem}_{iso}.xlsx", engine='xlsxwriter')
        df = df.T
        df.to_excel(writer, sheet_name='trace', startrow=1)
        wb = writer.book
        ws = writer.sheets['trace']
        ws.write(0, 0, 'trace')
        ws.write(1, 0, 'ppm')
        for i, item in enumerate(header[0]):
            ws.write(0, i + 1, item)
        if os.path.exists(f'{loc}/cfg_{iso}.txt'):
            # Write reference ppm config
            message("Writing reference peak config.", vb)
            cfgsheet = wb.add_worksheet('cfg')
            with open(f'{loc}/cfg_{iso}.txt', 'r') as rf:
                for i, line in enumerate(rf):
                    l = line.strip('\n').split('\t')
                    cfgsheet.write(i, 0, l[0])
                    cfgsheet.write(i, 1, l[1])
        else:
            message("Config file not found.", vb)
        traj.to_excel(writer, sheet_name='area', index=False)
        writer.save()
        message(f"Completed successfully. Output: {loc}/{stem}_{iso}.xlsx", vb)

def detect_spectra(iso='13C', init=11):
    """Returns a list of spectra IDs for the given isotope."""
    isotopes = iso.split("_")
    iso_spectra = []
    dirs = glob.glob('*/')
    st = [str(y) for y in sorted([int(x[:-1]) for x in dirs if str.isdigit(x[:-1])])]
    for fn in st:
        if int(fn) >= int(init):
            try:
                dic, data = ng.bruker.read(fn)
                # for i, label in enumerate(isotopes):
                #     match = True
                #     if dic['acqus'][f"NUC{i+1}"] != label:
                #         match = False
                # if match:
                #     iso_spectra.append(fn)
                udic = ng.bruker.guess_udic(dic, data)
                if (udic[0]['label'] == iso \
                        and udic[0]['encoding'] == 'direct' \
                        and dic['acqus']['PULPROG'] != 'zgpr_llc'):
                    iso_spectra.append(fn)
            except OSError:
                pass
    return iso_spectra

def main(iso, init):
    """Call stack processing function."""
    print("Welcome to NMR processing pipeline.\n")
    if iso not in ("1H", "13C", "1H_13C"):
        print("Isotope not recognized.")
        return
    loc = os.getcwd()
    # indices = detect_spectra(iso=iso, init=init)
    indices = [str(3*i + 1) for i in range(4, 34)]
    # indices = ['25']
    # indices = ['91'] #, '14', '15']
    # print("Beginning calibration with following parameters:")
    print("\t- Directory: " + loc)
    print(f"\t- {iso} spectrum IDs: " + ', '.join(indices))
    print("\t- Initial timepoint reference: " + init)
    process_stack(loc, indices, initial=init, iso=iso, vb=True)

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
    elif args[0] in ['13C', '1H', '31P', '1H_13C']:
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
