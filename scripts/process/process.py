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
import skimage

SCDIR = os.path.dirname(__file__)   # location of script
CURVE_FIT_RESOLUTION = 2**19

# ver = ''
# ver = '_2'
ver = '_3'
# ver = '_4'

ISOTOPE_PARAMS = {
    '1H': {
        'fit_ppm_delta': 0.1,
        'fit_fwhm_max': 0.2,
        'fit_fwhm_min': 0.001,
        'plot_bounds': (12., 0.),
    },
    '13C': {
        'fit_ppm_delta': 5,
        'fit_fwhm_max': 0.5,
        'fit_fwhm_min': 0.01,
        'plot_bounds': (200., 0.),
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
    if dic['FDF2LABEL'] == "1H":
        data = ng.process.proc_bl.baseline_corrector(data, wd=25)

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

        Parameters:
        reference_shifts -- set of reference ppm shifts corresponding to the
            overlapping peaks
        assignments -- dict mapping compounds to assigned subpeak indices

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

def vp_fit(dic, data, isotope, peaks, peak_width, cIDs=None):
    """Perform Voigt profile (Gaussian-Lorentzian convolution) fit for detected peaks.
    
    Parameters:
    dic -- NMR spectrum metadata to create UC object.
    data -- FT-NMR spectrum for fitting.
    isotope -- Primary isotope of experiment.
    peaks -- List of locations (data array indices) where peaks should be fit.
    peak_width -- Default width(s) of peaks. If a number, will be applied to all peaks.
                  If a list, must be same length as peaks and each width applies to
                  corresponding peak.
    cIDs -- optional, cluster IDs to group peaks together for fitting. Default: all fit
            together.

    Returns optimal parameters of fit peaks.
    """
    uc = ng.pipe.make_uc(dic, data, dim=0) # Create unit conversion object

    # Retreive reasonable parameter bounds for isotope
    ppm_delta = ISOTOPE_PARAMS[isotope]['fit_ppm_delta']
    fwhm_min = ISOTOPE_PARAMS[isotope]['fit_fwhm_min']
    fwhm_max = ISOTOPE_PARAMS[isotope]['fit_fwhm_max']
    amps = np.array([data.real[a] for a in peaks])
    amp_bounds = [(0.2*data.real[a], 1.5*data.real[a]) for a in peaks]

    # Handle peak widths
    if isinstance(peak_width, collections.abc.Sequence):
        peak_widths = peak_width
    else:
        peak_widths = [(peak_width,) for i in range(len(peaks))]

    # Handle cluster IDs
    if cIDs is None:
        cIDs = [1 for i in range(len(peaks))]

    # Define start parameters and bounds for each peak
    vf_params = [((a[0], 0.5*b[0], 0.5*b[0]),) for a, b in zip(peaks, peak_widths)]
    vf_bounds = []
    for par in vf_params:
        vf_bounds.append(((
            (par[0][0] - ppm_delta, par[0][0] + ppm_delta),    # shift bounds
            (ppm_ivl(uc, fwhm_min), ppm_ivl(uc, fwhm_max)), # gauss bounds
            (ppm_ivl(uc, fwhm_min), ppm_ivl(uc, fwhm_max)), # lorentz bounds
        ),))

    # Perform curve fit and return result
    params, amps, p_err, a_err, found = ng.analysis.linesh.fit_spectrum(
        data, ['v'], vf_params, amps, vf_bounds, amp_bounds, peaks, cIDs,
        int(ppm_ivl(uc, ppm_delta)), True, verb=False
    )
    return params, amps

def sg_peakpick(data, w_len, p_order, r=2):
    """Detect peaks using a second-derivative Savitzky-Golay filter.
    
    Parameters:
    data -- FT-NMR spectrum for peak-picking.
    w_len -- SG filter window size, must be odd. Larger window size results in greater
             smoothing.
    p_order -- order of SG polynomial.
    r -- number of standard deviations for peak-finding minimum prominence.
    
    Return array containing indices of found peaks.
    """
    fd2 = -1*scipy.signal.savgol_filter(data, w_len, p_order, deriv=2)
    err = r*np.std(fd2)
    pks, prp = scipy.signal.find_peaks(fd2, distance=3, width=(3, 75), prominence=err)
    x = np.array(range(len(data), 0, -1))
    pklabels = list(range(pks.shape[0]))
    plot_with_labels(x, fd2, pks, pklabels, ppm_bounds=(max(x), min(x)))
    return pks

def peak_pick(data, isotope='1H'):
    """Use Savitzky-Golay method to detect both fine and coarse peaks.
    
    TODO: tune for 13C
    
    Parameters:
    data -- vector of data for peak-picking.
    isotope -- isotope to determine SG filter parameters (currently unused,
               tuned to 1H.

    Returns list containing positions of peaks detected.
    """
    coarse_peaks = list(sg_peakpick(data, 31, 3))
    fine_peaks = list(sg_peakpick(data, 11, 2))
    peaks = [(pk,) for pk in sorted(fine_peaks + coarse_peaks)]
    return peaks

def peak_fit(dic, data, history, r=6, sep=0.005, plot=True, vb=False):
    """Peak-pick NMR spectrum.
    Includes peak-picking, peak assignment, and integration functionality. The
    reference peaks to be fit should be specified in cfg_{isotope}.txt and
    provided in the directory of the run to be processed.

    Parameters:
    ppm -- chemical shift axis of the data, in ppm
    data -- signal axis of the data
    history -- history of peak assignments in deconvolutions
    r -- RMSE noise threshold for peak detection (default: 6)
    plot -- True (default) to plot fit peaks with cluster labels over the data.
    vb -- True for verbose mode.

    Returns a dictionary mapping compounds to their peak areas, and a curve
    summing all of the fit peaks.
    """
    # ---------------------------- Utilities ----------------------------
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    data = data.real
    isotope = dic["FDF2LABEL"]
    plot_bounds = ISOTOPE_PARAMS[isotope]['plot_bounds']
    if isotope == '1H':
        data = -1*data

    def p2f_interval(xx): return ppm_ivl(uc, xx) # convert PPM scale to index
    def p2i_interval(xx): return int(p2f_interval(xx)) # convert PPM scale to index # convert PPM scale to index
    def p2i(xx): return uc.f(xx, unit='ppm')
    def i2p(yy): return uc.ppm(yy) # convert index to PPM shift
    def i2i(xx): return xx

    def message(label, *args):
        """Messages to print in verbose mode."""
        if not vb:
            return
        if label == 'closest_peak_message':
            text = "No reference peaks found within range of {0} ppm. The " \
                   "closest peak will be assigned:\n{1}"
        elif label == 'collision_message':
            text = "Attempted to assign the following reference peaks to " \
                   "cluster {0}, but a collision was detected.\n{1}"
        elif label == 'collision_reassign_message':
            text = "Cluster {0} is assigned to more than one peak, so its " \
                   "assignment at shift {1} will be overwritten by new " \
                   "cluster {2}."
        elif label == 'collision_abort_message':
            text = "Cluster {0} will be orphaned if the assignment for shift " \
                   "{1} is overwritten, so it will not be reassigned to " \
                   "cluster {2}."
        elif label == 'orphaned_message':
            text = "Cluster {0} is orphaned."
        print(text.format(*args))

    # ---------------------------- PEAK-PICKING -----------------------------
    if isotope == '1H':
        pthres = 1000 # 300
        shifts = peak_pick(data)
        cIDs = ng.analysis.peakpick.clusters(data, shifts, pthres, None, None, None, 
                                             p2i_interval(0.3))
        params, amps = vp_fit(dic, data, isotope, shifts, p2i_interval(0.1), cIDs=cIDs)
    else:
        pthres = r*rmse(data[(ppm <= 160) & (ppm >= 130)])
        shifts, cIDs, params, amps = ng.analysis.peakpick.pick( # Find peaks
            data, pthres=pthres, msep=(p2i_interval(sep),),
            algorithm='thres-fast', est_params=True, lineshapes=['g'],
            cluster=True, c_ndil=p2i_interval(0.3), table=False)

    amps = np.array([data.real[a] for a in shifts])
    if len(shifts) == 0:
        print("No peaks detected.")
        return {}, []
    shifts = np.array(shifts)

    if plot:
        plot_with_labels(ppm, data, shifts, cIDs, ppm_bounds=plot_bounds)
        plt.show()

    # -------------------- Filter noise artifacts on peaks -------------------
    def prominent(subpeaks, subamps, sep=0.1, r=0.1): # r=p2p(0.1) # sep=0.1, r=0.25
        """Finds peaks likely generated by noise distortions.
        These artifacts appear on some strong NMR signals around the base. This
        method subclusters close peaks separated by less than <sep> and filters
        out those less than <r> times the maximum subcluster amplitude.

        Parameters:
        subpeaks -- array of chemical shifts of peaks to be examined
        subamps -- array of amplitudes of peaks to be examined
        sep -- max distance for subclustering in ppm (default: 0.1)
        r -- threshold of peak height ratio for exclusion (default: 0.25)

        Returns a list of booleans representing "prominent" peaks to include
        from the cluster.
        """
        distances = scipy.spatial.distance.cdist(subpeaks, subpeaks)
        nclust, sub_cIDs = scipy.sparse.csgraph.connected_components(
            (distances <= sep).astype(int), 
            directed=False
        )
        include = []
        for i in set(sub_cIDs):
            ampset = subamps[sub_cIDs == i]
            include.extend(abs(ampset) >= r*abs(ampset).max(axis=0))
        return include

    def separated(data, subpeaks, subamps, sep=0.1, r=0.5, include=None):
        """Finds peaks likely generated by noise distortions.
        These artifacts appear on some strong NMR signals around the base. This
        method subclusters close peaks separated by less than <sep> and filters
        out those that are not separated by valleys of amplitude less than <r>
        times the maximum pair amplitude.

        Parameters:
        data -- array of spectrum data
        subpeaks -- array of chemical shifts of peaks to be examined
        subamps -- array of amplitudes of peaks to be examined
        sep -- max distance for subclustering in ppm (default: 0.1)
        r -- threshold of peak height ratio for exclusion (default: 0.25)

        Returns a list of booleans representing "prominent" peaks to include
        from the cluster.
        """
        if include is None:
            include = [True]*len(subpeaks)
        baseline = []
        while baseline != include:
            baseline = copy.deepcopy(include)
            idx_range = [i for i, v in enumerate(include) if v == True]
            for k, ii in enumerate(idx_range[:-1]):
                jj = idx_range[k + 1]
                p = subpeaks[ii]
                q = subpeaks[jj]
                a = subamps[ii]
                b = subamps[jj]
                if abs(q - p) < p2i_interval(sep):
                    idxs = (int(p), int(q))
                    if min(data[min(*idxs):max(*idxs)]) > r*min(a, b):
                        # absorb lower subpeak
                        if float(a) < float(b):
                            include[ii] = False
                        else:
                            include[jj] = False
            if not any(include):
                include[np.argmax(subamps)] = True
        return include

    mask = []
    for cid in np.unique(cIDs):
        include = cIDs == cid
        if isotope == '1H':
            submask = prominent(shifts[include], amps[include], sep=0.05, r=0.5)
            mask.extend(separated(data, shifts[include], amps[include], sep=0.05, r=0.5, include=submask))
        else:
            mask.extend(prominent(shifts[include], amps[include]))
    if isotope == '1H':
        mask = [True for a in shifts]
        mask = [False if i2p(shifts[i]) > 4.4 else e for i, e in enumerate(mask)]
    mask = np.array(mask)
    shifts = np.array(shifts)[mask]
    cIDs = np.array(cIDs)[mask]
    params = np.array(params)[mask]
    amps = np.array(amps)[mask]
    if len(shifts) == 0:
        plot_datasets(ppm, [data], ppm_bounds=plot_bounds)
        return {}, []
    elif plot:
        plot_with_labels(ppm, data, shifts, cIDs, ppm_bounds=plot_bounds)
    # Collapse clusters
    clust_map = {c: i for i, c in enumerate(np.unique(cIDs))}
    cIDs = np.array([clust_map[c] for c in cIDs])

    # Fit Voigt curves at found peaks
    # Guess lineshape parameters and bounds
    amp_bounds = [(0.8*data.real[a], 1.2*data.real[a]) for a in shifts]
    vf_params = [((a[0], 0.5*b[0], 0.5*b[0]),) for a, b in zip(shifts, params)]
    vf_bounds = []
    for par in vf_params:
        vf_bounds.append(((
            (par[0][0]-i2i(5), par[0][0]+i2i(5)),    # shift bounds
            (p2f_interval(0.01), p2f_interval(0.5)), # gauss bounds
            (p2f_interval(0.01), p2f_interval(0.5)), # lorentz bounds
        ),))
    params, amps, p_err, a_err, found = ng.analysis.linesh.fit_spectrum(
        data, ['v'], vf_params, amps, vf_bounds, amp_bounds, shifts, cIDs,
        p2i_interval(0.5), True, verb=False
    )
    mask = ~(np.isnan(amps) \
            | np.array([any(np.isnan(pset[0])) for pset in params]) \
            | (np.array([sh[0] for sh in shifts]) <= p2i(170)))
    shifts_tup = np.array(shifts)[mask]
    amp_bounds = np.array(amp_bounds)[mask]
    vf_params = np.array(vf_params)[mask]
    vf_bounds = np.array(vf_bounds)[mask]
    shifts = np.array([sh[0] for sh in shifts])[mask]
    params = np.array(params)[mask]
    cIDs = np.array(cIDs, dtype='l')[mask]
    amps = np.array(amps)[mask]
    if plot:
        """Plot found peaks."""
        plot_with_labels(ppm, data.real, shifts, cIDs, ppm_bounds=plot_bounds)

    # ---------------- Assign clusters to reference peaks --------------------
    refsh = pd.read_csv(f"cfg_{isotope}.txt", sep='\t', names=['Shift', 'Compound']) # Shift | Compound
    refsh['Assigned'] = np.nan # Create series for assigned clusters

    for cid in np.unique(cIDs):
        # Assign any reference peaks within the cluster bounds to this cluster
        subpeaks = np.array([i2p(sh) for sh in shifts[cIDs == cid]])
        mask = (refsh['Shift'] > min(subpeaks) - 0.45) \
                & (refsh['Shift'] < max(subpeaks) + 0.45) ## IMPORTANT: ADD CASE FOR 1H
        refpks = refsh.loc[mask, 'Shift'].tolist()

        # If no refpks w/in range, assign the nearest refpk w/in 5ppm to this cluster
        if len(refpks) == 0:
            cshift = i2p(subpeaks.mean())
            mask = (refsh['Shift']-cshift).abs().argsort()[0]
            nearest = refsh.iloc[mask, refsh.columns.get_loc('Shift')]
            if abs(nearest - cshift) < 5:
                refpks = [nearest]
                message('closest_peak_message', cshift, refpks)

        refrows = refsh['Shift'].isin(refpks)
        # Handle collisions:
        if any(~np.isnan(refsh.loc[refrows, 'Assigned'])):
            to_assign = refsh.loc[refrows, :]
            message('collision_message', cid, to_assign)
            colliding = to_assign.loc[~np.isnan(to_assign['Assigned']), :]
            for pk, row in colliding.iterrows():
                cl = colliding.at[pk, 'Assigned']
                if refsh.loc[refsh['Assigned'] == cl, :].shape[0] > 1:
                    message('collision_reassign_message', cl, pk, cid)
                else:
                    message('collision_abort_message', cl, pk, cid)
                    refpks.remove(colliding.at[pk, 'Shift'])

        if len(refpks) == 0:
            message('orphaned_message', cid)
        else:
            refsh.loc[refsh['Shift'].isin(refpks), 'Assigned'] = cid

    # ---------------------- Assign peaks to compounds -----------------------
    assignments = collections.defaultdict(list)
    for cid in np.unique(cIDs):
        # Assign cluster subpeaks to compounds
        include = cIDs == cid
        subamps = np.array(amps[include])
        subpeaks = np.array(shifts[include])
        assigned_ref = refsh[refsh['Assigned'] == cid] # reference shifts assigned to cluster
        cpds = np.unique(assigned_ref['Compound'])
        refpks = set(assigned_ref['Shift'].tolist())
        # Assign all subpeaks to one compound if it owns the cluster
        if cpds.shape[0] == 1:
            assignments[cpds[0]].extend(subpeaks)
        elif cpds.shape[0] > 1:
            # Assign using deconvolution history
            stored_pattern = False
            if history is not None:
                if history.has_pattern(refpks, len(subpeaks)):
                    stored_pattern = True
                    pattern = history.get_pattern(refpks, len(subpeaks))
                    for cpd, pks in pattern.get_assignments():
                        assignments[cpd].extend(subpeaks[pks])
            if not stored_pattern:
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
                                i2p(max(subpeaks))-0.45) ## IMPORTANT: ADD CASE FOR 1H
                )
                for cpd in cpds:
                    while True:
                        raw_text = input(f"Subpeak contributions of {cpd}: ")
                        try:
                            indices = [int(i) for i in raw_text.strip().split()]
                            break
                        except ValueError:
                            print("Invalid assignments.")
                    splitting_assignments[cpd] = indices
                    if len(indices) > 0:
                        assignments[cpd].extend(subpeaks[indices])
                if history is not None:
                    history.add_pattern(refpks, splitting_assignments)

    all_cpds = np.array([cpd for cpd, apks in assignments.items() if len(apks) > 0])

    # --------------------- CURVE-FIT PEAKS IN WINDOWS -----------------------
    def curve_fit(dic, data, shifts, vf_params, vf_bounds, amps, amp_bounds,
                  cIDs):
        """Perform nmrglue Voigt curve fit on data with parameters."""
        params, amps, found = ng.analysis.linesh.fit_spectrum(
            data, ['v'], vf_params, amps, vf_bounds, amp_bounds, shifts, cIDs,
            p2i_interval(0.5), False, verb=False
        )
        return params, amps, found
    # Extract islands with peaks separated by <= 10ppm, including 8ppm buffer
    curve_params = [shifts_tup, vf_params, vf_bounds, amps, amp_bounds, cIDs]
    if isotope == '1H':
        fit_params, fit_amps, found = curve_fit(dic, data, *curve_params)
    else:
        fit_params = []
        fit_amps = []
        p = 0
        q = 0
        for i, shift in enumerate(shifts):
            last = len(shifts) - 1
            if i == last:
                q = i
            if (abs(shift - shifts[q]) > p2i_interval(10)) or (i == last):
                # extract region from p-8 to q+8
                lb = shifts[p] - p2f_interval(8)
                ub = shifts[q] + p2f_interval(8)
                ndic, ndata = ng.process.pipe_proc.ext(dic, data, x1=lb, xn=ub)

                # Shift peak parameters and perform curve fit
                args = [a[p:q+1] for a in curve_params]
                args[0] = shift_params(lb, args[0])[0]
                args[1] = shift_params(lb, *args[1])
                args[2] = shift_params(lb, *args[2])
                nparams, namps, found = curve_fit(ndic, ndata, *args)
                if False: # plot fit peaks with cluster labels over the data.
                    nppm = ng.pipe.make_uc(ndic, ndata, dim=0).ppm_scale()
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(nppm, ndata, lw=0.5)
                    sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
                        ndata.shape, ['v'], nparams, namps)
                    ax.plot(nppm, sim, lw=0.5)
                    plt.show()
                nparams = shift_params(-1*lb, *np.array(nparams))
                fit_params.extend(nparams)
                fit_amps.extend(list(namps))
                p = i
            q = i
    fit_params = np.array(fit_params)
    fit_amps = np.array(fit_amps)
    resolution = CURVE_FIT_RESOLUTION
    scale = resolution/data.shape[0]
    sppm = np.linspace(ppm[0], ppm[-1], num=resolution)
    if plot: # plot fit peaks with cluster labels over the data.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5)
        for cpd in all_cpds:
            color = next(ax._get_lines.prop_cycler)['color']
            include = np.isin(shifts, assignments[cpd])
            idxs = shifts[include]
            sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
                (resolution,), ['v'], scale*fit_params[include], fit_amps[include])
            ax.plot(sppm, sim, lw=0.5, color=color)
            ax.plot(ppm[idxs], fit_amps[include], ls='', marker='*', color=color)
            for xx, yy in zip(ppm[idxs], fit_amps[include]):
                ax.annotate(cpd, (xx,yy), textcoords='offset points',
                    xytext=(0,10), ha='center')
        ax.set_xlim(*plot_bounds)
        # ax.set_xlim(2.5, 0.5)
        # ax.set_ylim(-250, 2500)
        plt.show()

    # --------------------------- INTEGRATE PEAKS ----------------------------
    # for each compound, integrate the simulated signal of just the peaks
    # belonging to that compound
    signals = {}
    curve = np.zeros((resolution,))
    for cpd in all_cpds:
        include = np.isin(shifts, assignments[cpd])
        sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
            (resolution,), ['v'], scale*fit_params[include], fit_amps[include])
        areacurve = -1*np.trapz(sim, x=sppm)
        signals[cpd] = areacurve
        curve += sim

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

def process_trace(loc, item, cf, ps, history, overwrite=True, man_ps=False):
    """Process a single 1D 13C-NMR trace.

    Parameters:
    loc -- path to folder containing traces
    item -- ID of trace to process
    cf -- correction factor for ppm shift
    ps -- list with the default p0, p1 values for phase correction
    history -- history of peak assignments in deconvolutions
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
        signals, curve = peak_fit(dic, data, history, r=4, plot=False, vb=True)
        # signals = ridge_trace(dic, data, plot=True)
        # curve = data.real
    else:
        signals, curve = peak_fit(dic, data, history, r=12, sep=30/2*0.005, plot=True, vb=True) # r=4 # sep=30/2*0.001
        # signals = ridge_trace(dic, data, plot=True)
        curve = data.real

    return ppm, data.real, timestamp, signals, n_scans, curve

def message(message, verbose):
    """If in verbose mode, print status message"""
    if verbose:
        print(message)

def process_stack(loc, ids, initial=None, iso='13C', vb=True, dry_run=True):
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
    # history = DeconvolutionHandler()
    history = None
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
        ppm, trace, ts, signals, n_scans, curve = process_trace(
            loc, item, cf, phases, history
        )
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
    # indices = ['94'] #, '14', '15']
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
