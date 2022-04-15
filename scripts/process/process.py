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
import contextlib
import copy
import datetime
import glob
import io
import os
import subprocess
import sys
import time

import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np
import pandas as pd

SCDIR = os.path.dirname(__file__)   # location of script

ver = '_2' # '_2'  # use '' or '_2'

def plot_datasets(xx, yys, ppm_bounds=(200., 0.)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for yy in yys:
        ax.plot(xx, yy, lw=0.5)
    ax.set_xlim(*ppm_bounds)
    plt.show()

def plot_with_labels(xx, yy, labeli, labels, ppm_bounds=(200., 0.)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xx, yy, lw=0.5)
    ax.plot(xx[labeli], yy[labeli], ls='', marker='*')
    for ci, xi, yi in zip(labels, xx[labeli], yy[labeli]):
        ax.annotate(ci, (xi,yi), textcoords='offset points', xytext=(0,10),
            ha='center')
    ax.set_xlim(*ppm_bounds)
    plt.show()

def find_ps(dic, data, ppm_bounds=(200., 0.)):
    ## CALIBRATE PHASE CORRECTION IF NECESSARY ##
    p0 = 0.
    p1 = 0.
    while True:
        print(f"Current phase correction: p0={p0}, p1={p1}")
        dic1, data1 = ng.process.pipe_proc.ps(dic, data, p0=p0, p1=p1)
        uc = ng.pipe.make_uc(dic1, data1, dim=0)
        plot_datasets(uc.ppm_scale(), [data1], ppm_bounds=ppm_bounds)
        text = input("Enter new p0 and p1 separated by whitespace, or 'done': ")
        if text in ['done', 'exit', 'quit', 'q', '', "'done'"]:
            break
        vals = text.split()
        try:
            p0_try = float(vals[0])
            p1_try = float(vals[1])
            p0 = p0_try
            p1 = p1_try
        except ValueError:
            print("Invalid shift string.")
    return p0, p1

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

def calibrate_process(loc, item, acq1, isotope, verbose=True):
    """Get process parameters using NMRglue and write NMRpipe script

    Parameters:
    loc -- directory for the NMR run
    item  -- ID of spectrum (also folder name)
    acq1 -- Initial acquisition ID for timepoint anchor

    Returns: ppm correction to calibrate the x axis
    """
    if isotope == '13C':
        ppm_bounds = (200., 0)
    elif isotope == '1H':
        ppm_bounds = (12., 0.)

    runf = loc.split('/')[-1]
    ## GET INITIAL TIMESTAMP ##
    time0 = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{acq1}/pulseprogram"))
    message("Loaded initial timestamp.", verbose)

    ## GET BRUKER PARAMS (Number of Scans) ##
    dic, data = ng.bruker.read(f"{loc}/{item}")
    n_scans = dic['acqus']['NS']

    ## CONVERT, PROCESS, and LOAD SPECTRUM ##
    pipe_output = subprocess.run(
        ["csh", f"{SCDIR}/proc_{isotope}a.com", item],
        stdout=subprocess.PIPE)
    dic, data = ng.fileio.pipe.read(pipe_output.stdout)

    ## LOAD CALIBRATION PARAMETERS AND PHASE SHIFT
    calibrations = pd.read_csv(f"{SCDIR}/calibrations.txt", sep='\t', dtype='str')
    rundate = runf.split('_')[0]
    calibr_entry = calibrations.loc[(calibrations["run"] == rundate) \
                                    & (calibrations["isotope"] == isotope)]
    run_calibrated = calibr_entry.shape[0] >= 1
    if not run_calibrated:
        p0, p1 = find_ps(dic, data, ppm_bounds=ppm_bounds)
    else:
        if calibr_entry.shape[0] > 1:
            print("More than one calibration entry matches this run; " \
                  "taking just the first.")
        entry_index = calibr_entry.index[0]
        ref_shift = float(calibr_entry.at[entry_index, 'shift_reference'])
        actual_shift = float(calibr_entry.at[entry_index, 'shift_actual'])
        p0 = float(calibr_entry.at[entry_index, 'p0'])
        p1 = float(calibr_entry.at[entry_index, 'p1'])
    dic, data = ng.process.pipe_proc.ps(dic, data, p0=p0, p1=p1)
    dic, data = ng.process.pipe_proc.mult(dic, data, r=n_scans, inv=True)
    data = ng.process.proc_bl.baseline_corrector(data, wd=20)

    ## CALIBRATE PPM SHIFT IF NECESSARY ##
    if not run_calibrated:
        ref_shift = 72.405
        print("Find the actual ppm shift of the glucose peak with" \
            f"reference shift {ref_shift}.")
        print("Note the x-coordinate, then close the window.")
        uc = ng.pipe.make_uc(dic, data, dim=0)
        plot_datasets(uc.ppm_scale(), [data.real], ppm_bounds=ppm_bounds)
        calib = input("Type the ppm shift of the peak you wish to calibrate: ")
        actual_shift = float(calib)

    calibration_shift = float(ref_shift) - float(actual_shift)
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    return calibration_shift, time0, [p0, p1], ppm

class DeconvolutionHandler:
    """Handles the history of peak splitting patterns."""
    def __init__(self):
        self._patterns = set([])

    def has_pattern(self, reference_shifts, subpeaks):
        """Determine whether a matching pattern exists in the history."""
        for pattern in self._patterns:
            if pattern.is_equivalent(reference_shifts, subpeaks):
                return True
        return False

    def add_pattern(self, reference_shifts, assignments):
        """Construct and a new deconvolution pattern or return an existing one.

        Parameters
        ----------
        reference_shifts: set of reference ppm shifts corresponding to the
            overlapping peaks
        assignments : dict mapping compounds to assigned subpeak indices

        Returns: the constructed pattern.
        """
        new_deconvolution = Deconvolution(reference_shifts, assignments)
        new_pattern = (new_deconvolution.num_subpeaks(),
                       new_deconvolution.get_reference_shifts())
        if self.has_pattern(*new_pattern):
            return self.get_pattern(*new_pattern)
        else:
            self._patterns.add(new_deconvolution)
            new_deconvolution.assign_handler(self)
            return new_deconvolution

    def get_pattern(self, reference_shifts, n_subpeaks):
        """Retrieve an existing deconvolution pattern."""
        for pattern in self._patterns:
            if pattern.is_equivalent(reference_shifts, n_subpeaks):
                return pattern
        return None

class Deconvolution:
    """Contains information about a peak splitting pattern."""
    def __init__(self, reference_shifts, assignments):
        """Construct a new deconvolution pattern with peak assignments."""
        self._reference_shifts = set(reference_shifts)
        self._assignments = dict(assignments)
        self._compounds = set([compound for compound in self._assignments])
        self._subpeaks = set([])
        for subpeak_set in self._assignments.values():
            self._subpeaks |= set(subpeak_set)
        self._handler = None

    def assign_handler(self, handler):
        self._handler = handler

    def num_subpeaks(self):
        return len(self._subpeaks)

    def get_reference_shifts(self):
        return self._reference_shifts

    def is_equivalent(self, reference_shifts, n_subpeaks):
        """Determine whether this pattern is equivalent to a query.
        An equivalent pattern has the same assigned reference peaks and the same
        number of subpeaks."""
        return (self._reference_shifts == reference_shifts) and (self.num_subpeaks() == n_subpeaks)

    def get_assignments(self):
        for compound, subpeaks in self._assignments.items():
            yield compound, subpeaks

def rmse(vector):
    """Calculate the Root Mean Squared Error (RMSE) of a 1D-array."""
    return np.sqrt(((vector - vector.mean()) ** 2).mean())

def peak_fit(dic, data, history, plot=False, vb=False):
    """Peak-pick NMR spectrum.
    Includes peak-picking, peak assignment, and integration functionality.

    Parameters:
    ppm -- chemical shift axis of the data, in ppm
    data -- signal axis of the data
    history -- history of peak assignments in deconvolutions
    plot -- True (default) to plot fit peaks with cluster labels over
        the data.
    vb -- True for verbose mode.
    """
    # Utilities
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    data = data.real
    def p2i_interval(xx): return int(abs(uc.f(0, unit='ppm') - uc.f(xx, unit='ppm'))) # convert PPM scale to index
    def p2i(xx): return uc.i(xx, unit='ppm')
    def i2p(yy): return uc.ppm(yy) # convert index to PPM shift
    def message(message):
        if vb:
            print(message)
    # pthres = 2*max(data[(ppm <= 160) & (ppm >= 130)])
    pthres = 10*rmse(data[(ppm <= 160) & (ppm >= 130)])
    shifts, cIDs, params, amps = ng.analysis.peakpick.pick( # Find peaks
        data, pthres=pthres, msep=(p2i_interval(0.08),), algorithm='thres',
        est_params=True, lineshapes=['g'], cluster=True,
        c_ndil=p2i_interval(0.3), table=False)
    shifts = np.array([sh[0] for sh in shifts])
    cIDs = np.array(cIDs, dtype='l')
    amps = np.array([max(data[(ppm < i2p(s) + 0.45) & (ppm > i2p(s) - 0.45)]) \
                     for s in shifts])
    if plot:
        """Plot found peaks."""
        plot_with_labels(ppm, data.real, shifts, cIDs)

    # Assign peaks
    refsh = pd.read_csv("cfg.txt", sep='\t', names=['Shift', 'Compound']) # Shift | Compound
    refsh['Assigned'] = np.nan # Create series for assigned clusters

    for cid in np.unique(cIDs):
        # Assign any reference peaks within the cluster bounds to this cluster
        subpeaks = np.array([i2p(sh) for sh in shifts[cIDs == cid]])
        mask = (refsh['Shift'] > min(subpeaks) - 0.45) \
                & (refsh['Shift'] < max(subpeaks) + 0.45)
        refpks = refsh.loc[mask, 'Shift'].tolist()

        # If no refpks w/in range, assign the nearest refpk to this cluster
        if len(refpks) == 0:
            mask = (refsh['Shift']-i2p(subpeaks.mean())).abs().argsort()[0]
            refpks = [refsh.iloc[mask, refsh.columns.get_loc('Shift')]]
            message("No reference peaks found within range of " \
                    f"{i2p(subpeaks.mean())} ppm.")
            message("The closest peak will be assigned:")
            message(refpks)
        refrows = refsh['Shift'].isin(refpks)
        # Handle collisions:
        if any(~np.isnan(refsh.loc[refrows, 'Assigned'])):
            message("Attempted to assign the following reference peaks to " \
                    f"cluster {cid}, but a collision was detected.")
            to_assign = refsh.loc[refrows, :]
            message(to_assign)
            colliding = to_assign.loc[~np.isnan(to_assign['Assigned']), :]
            for pk, row in colliding.iterrows():
                cl = colliding.at[pk, 'Assigned']
                if refsh.loc[refsh['Assigned'] == cl, :].shape[0] > 1:
                    message(f"Cluster {cl} is assigned to more than one peak," \
                            f" so its assignment at shift {pk} will be" \
                            f" overwritten by new cluster {cid}.")
                else:
                    message(f"Cluster {cl} will be orphaned if the assignment" \
                            f" for shift {pk} is overwritten, so it will not" \
                            f" be reassigned to cluster {cid}.")
                    refpks.remove(colliding.at[pk, 'Shift'])

        if len(refpks) == 0:
            message(f"Cluster {cid} is orphaned.")
        else:
            refsh.loc[refsh['Shift'].isin(refpks), 'Assigned'] = cid

    # Correct reference shifts using assigned peaks
    assignment_map = {}
    def assign_subpeaks(ref, cpd, subpeaks, map=assignment_map):
        reference_peak = ref.loc[ref['Compound'] == cpd, 'Shift']
        refpk_index = map[cpd].index(p2i(reference_peak))
        map[cpd][refpk_index] = subpeaks.mean()

    for cpd in np.unique(refsh['Compound']):
        # Populate peak assignment map with reference shifts
        reference_peaks = refsh.loc[refsh['Compound'] == cpd, 'Shift']
        assignment_map[cpd] = [p2i(shift) for shift in reference_peaks]
    for cid in np.unique(cIDs):
        # Assign cluster subpeaks to compounds
        subpeaks = np.array(shifts[cIDs == cid])
        assigned_ref = refsh[refsh['Assigned'] == cid] # reference shifts assigned to cluster
        cpds = np.unique(assigned_ref['Compound'])
        # Assign all subpeaks to one compound if it owns the cluster
        if cpds.shape[0] == 1:
            assign_subpeaks(assigned_ref, cpds[0], subpeaks)
        else:
            refpks = set(assigned_ref['Shift'].tolist())
            # Assign using deconvolution history
            if history.has_pattern(refpks, len(subpeaks)):
                pattern = history.get_pattern(refpks, len(subpeaks))
                for cpd, pks in pattern.get_assignments():
                    assign_subpeaks(assigned_ref, cpd, subpeaks[pks])
            else:
                if len(refpks) > 0:
                    splitting_assignments = {}
                print(f"Cluster {cid} has {cpds.shape[0]} compounds.")
                print("In the plot, find the contributions from each " \
                      "compound using the reference shifts listed:")
                with pd.option_context(
                        'display.max_rows', None, 'display.max_columns', None):
                    print(assigned_ref)
                plot_with_labels(
                    ppm, data.real, subpeaks, range(len(subpeaks)),
                    ppm_bounds=(i2p(min(subpeaks))+0.45,
                                i2p(max(subpeaks))-0.45)
                )
                for cpd in cpds:
                    while True:
                        raw_text = input(f"Subpeak contributions of {cpd}: ")
                        try:
                            indices = [int(i) for i in raw_text.strip().split()]
                            break
                        except ValueError:
                            print("Invalid assignments.")
                    assign_subpeaks(assigned_ref, cpd, subpeaks[indices])
                    if len(indices) > 0:
                        splitting_assignments[cpd] = indices
                history.add_pattern(refpks, splitting_assignments)

    all_cpds = np.unique(refsh.loc[~np.isnan(refsh['Assigned']), 'Compound'])

    # mask = ~np.isnan(refsh['Assigned'])
    if plot: # plot fit peaks with cluster labels over the data.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5)
        for cpd in all_cpds:
            mask = np.full((data.shape[0],), False)
            color = next(ax._get_lines.prop_cycler)['color']
            idxs = [int(sh) for sh in assignment_map[cpd]]
            for sh in idxs:
                mask = mask | (ppm <= i2p(sh) + 0.45) & (ppm >= i2p(sh) - 0.45)
            ax.plot(ppm[mask], data[mask], lw=1, color=color)
            hts = np.array([max(data[(ppm < i2p(s) + 0.45) & (ppm>i2p(s)-0.45)]) for s in idxs])
            ax.plot(ppm[idxs], hts, ls='', marker='*', color=color)
            for xx, yy in zip(ppm[idxs], hts):
                ax.annotate(cpd, (xx,yy), textcoords='offset points',
                    xytext=(0,10), ha='center')
        ax.set_xlim(200, 0)
        plt.show()
    # Integrate Peaks - for each compound, integrate the simulated signal of
    # just the peaks belonging to that compound
    signals = {}
    for cpd in all_cpds:
        mask = np.full((data.shape[0],), False)
        idxs = assignment_map[cpd]
        for sh in idxs:
            mask = mask | (ppm <= i2p(sh) + 0.45) & (ppm >= i2p(sh) - 0.45)
        cpdsignal = copy.deepcopy(data)
        cpdsignal[~mask] = 0.
        areashape = ng.analysis.integration.integrate(cpdsignal, uc, (0, 200))
        signals[cpd] = areashape[-1]
    return signals

def process_trace(loc, item, cf, phases, isotope, history):
    """Process a single 1D 13C-NMR trace.

    Parameters:
    loc -- path to folder containing traces
    item -- ID of trace to process
    cf -- correction factor for ppm shift
    phases --
    isotope --
    history -- history of peak assignments in deconvolutions

    Returns:
    ppm --  list of ppm values
    trace -- list of signal intensities
    tim --  timestamp of trace, calculated from acquisition data
    phases -- P0, P1 phase shift
    """
    ## GET TIMESTAMP ##
    start_time = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{item}/pulseprogram"))
    end_time = datetime.datetime.fromtimestamp(
        os.path.getmtime(f"{loc}/{item}/acqus"))
    timestamp = start_time + (end_time - start_time)/2

    ## GET BRUKER PARAMS (NS) ##
    bruker_dic, _ = ng.bruker.read(f"{loc}/{item}")
    n_scans = bruker_dic['acqus']['NS']

    ## CONVERT, PROCESS, and LOAD SPECTRUM ##
    pipe_output = subprocess.run(
        ["csh", f"{SCDIR}/proc_{isotope}a.com", item],
        stdout=subprocess.PIPE)
    dic, data = ng.fileio.pipe.read(pipe_output.stdout)

    ## PHASE SHIFT, BASELINE CORRECTION, NORMALIZATION, and CALIBRATION ##
    dic, data = ng.process.pipe_proc.ps(dic, data, p0=phases[0], p1=phases[1])
    # data = ng.process.proc_bl.baseline_corrector(data, wd=20)

    # data = ng.process.proc_autophase.autops(data, 'peak_minima', p0=phases[0], p1=phases[1], peak_width=150)
    dic, data = ng.process.pipe_proc.mult(dic, data, r=n_scans, inv=True)
    pipe_output = subprocess.run(
        ["csh", f"{SCDIR}/proc_{isotope}b.com"],
        input=to_stream(dic, data),
        stdout=subprocess.PIPE)
    dic, data = ng.fileio.pipe.read(pipe_output.stdout)
    data = ng.process.proc_bl.baseline_corrector(data, wd=60)
    ppm = ng.pipe.make_uc(dic, data, dim=0).ppm_scale()
    shift = cf*data.shape[0]/(max(ppm) - min(ppm))
    dic, data = ng.process.pipe_proc.ls(dic, data, shift, sw=False)

    # Pick, fit, and integrate peaks
    signals = peak_fit(dic, data, history, plot=True)
    ppm = ng.pipe.make_uc(dic, data, dim=0).ppm_scale()

    return ppm, data.real, timestamp, signals, n_scans

def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, isotope='13C', c_each=False, verbose=True):
    """Process entire NMR stack and write to xlsx file.

    Parameters:
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
    history = DeconvolutionHandler()

    ## calibrate and get correction factor ##
    message("Calibrate autophase and PPM shift...", verbose)
    cf, time0, phases, c_ppm = calibrate_process(loc, ids[0], initial, isotope, verbose)
    c_ppm = c_ppm[(c_ppm <= 200) & (c_ppm >= 0)]
    stack.append(c_ppm)
    ## process all traces ##
    message("Begin processing " + isotope + " traces.", verbose)
    for i, item in enumerate(ids):
        message(f"Processing file {item}...", verbose)
        if c_each: # Calibrate spectra individually if -ce flag
            cf, _, phases, _ = calibrate_process(loc, item, initial, isotope, verbose)

        # Process spec and chop signal vector
        ppm, trace, ts, signals, n_scans = process_trace(loc, item, cf, phases, isotope, history)
        start_index = next(i for i, shift in enumerate(ppm) if shift < 200)
        idx = (start_index, start_index+len(stack[0]))
        vector = trace[idx[0]:idx[1]]

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
    traj.to_excel(writer, sheet_name='area', index=False)
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
                dic, data = ng.bruker.read(fn)
                udic = ng.bruker.guess_udic(dic, data)
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
