"""Old version of peak-fitting algorithm"""

def peak_fit(dic, data, history, plot=True):
    """Peak-pick NMR spectrum with Voigt lineshapes.
    Includes peak-picking, curve-fitting, and peak assignment functionality.

    Parameters:
    ppm -- chemical shift axis of the data, in ppm
    data -- signal axis of the data
    history -- history of peak assignments in deconvolutions
    plot -- True (default) to plot fit peaks with cluster labels over
        the data.
    """
    uc = ng.pipe.make_uc(dic, data, dim=0)
    ppm = uc.ppm_scale()
    data = data.real
    def ppc(xx): return int(abs(uc.f(0, unit='ppm') - uc.f(xx, unit='ppm'))) # convert PPM scale to index
    def ppsh(xx): return uc.i(xx, unit='ppm')
    def ixc(yy): return uc.ppm(yy) # convert index to PPM shift
    # pthres = 2*max(data[(ppm <= 160) & (ppm >= 130)])
    pthres = 10*rmse(data[(ppm <= 160) & (ppm >= 130)])
    shifts, cIDs, params, amps = ng.analysis.peakpick.pick( # Find peaks
        data, pthres=pthres, msep=(ppc(0.08),), algorithm='thres',
        est_params=True, lineshapes=['g'], cluster=True, c_ndil=ppc(0.3),
        table=False)

    # Fit Voigt curves at found peaks
    amps = [data.real[a] for a in shifts]
    amp_bounds = [(0.8*data.real[a], 1.2*data.real[a]) for a in shifts]
    vf_params = [((a[0], 0.5*b[0], 0.5*b[0]),) for a, b in zip(shifts, params)]
    vf_bounds = [(((par[0][0]-5, par[0][0]+5), (0.01, ppc(0.5)),
        (0.01, ppc(0.5))),) for par in vf_params]
    params, amps, p_err, a_err, found = ng.analysis.linesh.fit_spectrum(
        data, ['v'], vf_params, amps, vf_bounds, amp_bounds, shifts, cIDs,
        ppc(0.5), True, verb=False)
    # mask = ~(np.isnan(amps) | np.array([any(np.isnan(pset[0])) for pset in params]))
    shifts = np.array([sh[0] for sh in shifts])#[mask]
    # params = np.array(params)[mask]
    cIDs = np.array(cIDs, dtype='l')#[mask]
    # amps = np.array(amps)[mask]
    amps = np.array([max(data[(ppm<ixc(s)+0.45)&(ppm>ixc(s)-0.45)]) for s in shifts])

    if plot: # plot fit peaks with cluster labels over the data.
        # print(cIDs)
        # print(shifts)
        # print(params)
        # print(amps)
        # sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
        #     data.shape, ['v'], params, amps)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5)
        # print(shifts)
        # ax.plot(ppm, sim, lw=0.5)
        ax.plot(ppm[shifts], data.real[shifts], ls='', marker='*')
        for ci, xx, yy in zip(cIDs, ppm[shifts], data.real[shifts]):
            ax.annotate(ci, (xx,yy), textcoords='offset points', xytext=(0,10),
                ha='center')
        ax.set_xlim(200, 0)
        plt.show()
        plt.close('all')

    # Assign peaks
    # for each cluster, determine the reference peak(s) associated with that cluster
    # if the cluster has just one reference peak, assign all cluster subpeaks to that compound
    # if the cluster has more than one reference peak associated with it:
    #   - plot the simulated spectrum over the raw spectrum, with subpeak labels
    #   - get a list of subpeaks associated with each reference peak
    #   - assign the subpeaks accordingly
    #   - store the subpeak assignment for future deconvolutions
    refsh = pd.read_csv("cfg.txt", sep='\t', names=['Shift', 'Compound']) # Shift | Compound
    refsh['Assigned'] = np.full((refsh.shape[0],), np.nan)
    for cid in np.unique(cIDs):
        # assign any reference peaks within the cluster bounds to this cluster
        subpeaks = np.array([ixc(sh) for sh in shifts[cIDs == cid]])
        minsh = min(subpeaks) - 0.45
        maxsh = max(subpeaks) + 0.45
        mask = (refsh['Shift'] > minsh) & (refsh['Shift'] < maxsh)
        refpks = refsh.loc[mask, 'Shift'].tolist()
        # if no refpks w/in range, assign the nearest refpk to this cluster
        if len(refpks) == 0:
            mask = (refsh['Shift']-ixc(subpeaks.mean())).abs().argsort()[0]
            refpks = [refsh.iloc[mask, refsh.columns.get_loc('Shift')]]
            # print(f"No reference peaks found within range of {ixc(subpeaks.mean())} ppm.")
            # print("The closest peak will be assigned:")
            # print(refpks)
        refrows = refsh['Shift'].isin(refpks)
        if any(~np.isnan(refsh.loc[refrows, 'Assigned'])):
            print("Attempted to assign the following reference peaks to " \
                  f"cluster {cid}, but a collision was detected.")
            to_assign = refsh.loc[refrows, :]
            print(to_assign)
            colliding = to_assign.loc[~np.isnan(to_assign['Assigned']), :]
            for pk, row in colliding.iterrows():
                cl = colliding.at[pk, 'Assigned']
                if refsh.loc[refsh['Assigned'] == cl, :].shape[0] > 1:
                    print(f"Cluster {cl} is assigned to more than one peak, " \
                          f"so its assignment at shift {pk} will be " \
                          f"overwritten by new cluster {cid}.")
                else:
                    print(f"Cluster {cl} will be orphaned if the assignment " \
                          f"for shift {pk} is overwritten, so it will not " \
                          f"be reassigned to cluster {cid}.")
                    refpks.remove(colliding.at[pk, 'Shift'])
        if len(refpks) == 0:
            print(f"Cluster {cid} is orphaned.")
        else:
            refsh.loc[refsh['Shift'].isin(refpks), 'Assigned'] = cid

    assignments = collections.defaultdict(list)
    assignment_map = {}
    unique_cpds = np.unique(refsh['Compound'])
    for cpd in unique_cpds:
        reference_peaks = refsh.loc[refsh['Compound'] == cpd, 'Shift']
        print([shift for _, shift in reference_peaks.items()])
        assignment_map[cpd] = [ppsh(shift) for _, shift in reference_peaks.items()]
        print([ixc(shift) for shift in assignment_map[cpd]])
    for cid in np.unique(cIDs):
        subpeaks = np.array(shifts[cIDs == cid])
        assigned = refsh[refsh['Assigned'] == cid]
        cpds = np.unique(assigned['Compound'])
        if cpds.shape[0] == 1:
            assignments[cpds[0]].extend(subpeaks)
            reference_peak = assigned.iat[0, assigned.columns.get_loc('Shift')]
            refpk_index = assignment_map[cpds[0]].index(ppsh(reference_peak))
            print(f"{cpds}: {assignment_map[cpds[0]][refpk_index]} to {subpeaks.mean()}")
            assignment_map[cpds[0]][refpk_index] = subpeaks.mean()
        else:
            refpks = set(assigned['Shift'].tolist())
            if history.has_pattern(refpks, len(subpeaks)):
                pattern = history.get_pattern(refpks, len(subpeaks))
                for cpd, pks in pattern.get_assignments():
                    assignments[cpd].extend(subpeaks[pks])
                    reference_peak = assigned.loc[refsh['Compound'] == cpd, 'Shift']
                    refpk_index = assignment_map[cpd].index(ppsh(reference_peak))
                    assignment_map[cpd][refpk_index] = subpeaks[pks].mean()
            else:
                if len(refpks) > 0:
                    splitting_assignments = {}
                print(f"Cluster {cid} has {cpds.shape[0]} compounds.")
                print("In the plot, find the contributions from each " \
                      "compound using the reference shifts listed:")
                with pd.option_context(
                        'display.max_rows', None, 'display.max_columns', None):
                    print(assigned)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(ppm, data.real, lw=0.5)
                # sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
                #     data.shape, ['v'], params[cIDs == cid], amps[cIDs == cid])
                # ax.plot(ppm, sim, lw=0.5)
                x = ppm[subpeaks]
                y = data.real[subpeaks]
                ax.plot(x, y, ls='', marker='*')
                for i, coords in enumerate(zip(x, y)):
                    ax.annotate(i, coords, textcoords='offset points',
                        xytext=(0,10), ha='center')
                ax.set_xlim(ixc(min(subpeaks))+0.45, ixc(max(subpeaks))-0.45)
                plt.show()
                plt.close('all')
                for cpd in cpds:
                    while True:
                        raw_text = input(f"Subpeak contributions of {cpd}: ")
                        try:
                            indices = [int(i) for i in raw_text.strip().split()]
                            break
                        except ValueError:
                            print("Invalid assignments.")
                    assignments[cpd].extend(subpeaks[indices])
                    reference_peak = assigned.loc[refsh['Compound'] == cpd, 'Shift']
                    refpk_index = assignment_map[cpd].index(ppsh(reference_peak))
                    assignment_map[cpd][refpk_index] = subpeaks[indices].mean()
                    if len(indices) > 0:
                        splitting_assignments[cpd] = indices
                history.add_pattern(refpks, splitting_assignments)

    all_cpds = np.unique(refsh.loc[~np.isnan(refsh['Assigned']), 'Compound'])
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5)
        for cpd in all_cpds:
            ran = np.full((data.shape[0],), False)
            include = np.isin(shifts, assignments[cpd])
            idxs = shifts[include]
            for sh in idxs:
                ran = ran | (ppm <= ixc(sh) + 0.45) & (ppm >= ixc(sh) - 0.45)
            ax.plot(ppm[ran], data[ran], lw=1)

            ax.set_xlim(200, 0)
        plt.show()
        plt.close('all')

    # mask = ~np.isnan(refsh['Assigned'])
    if plot: # plot fit peaks with cluster labels over the data.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ppm, data.real, lw=0.5)
        for cpd in all_cpds:
            ran = np.full((data.shape[0],), False)
            color = next(ax._get_lines.prop_cycler)['color']
            include = np.isin(shifts, assignments[cpd])
            idxs = shifts[include]
            idxs = [int(sh) for sh in assignment_map[cpd]]
            for sh in idxs:
                ran = ran | (ppm <= ixc(sh) + 0.45) & (ppm >= ixc(sh) - 0.45)
            ax.plot(ppm[ran], data[ran], lw=1, color=color)
            # sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
            #     data.shape, ['v'], params[include], amps[include])
            # ax.plot(ppm, sim, lw=0.5, color=color)
            hts = np.array([max(data[(ppm<ixc(s)+0.45)&(ppm>ixc(s)-0.45)]) for s in idxs])
            ax.plot(ppm[idxs], hts, ls='', marker='*', color=color)
            for xx, yy in zip(ppm[idxs], hts):
                ax.annotate(cpd, (xx,yy), textcoords='offset points',
                    xytext=(0,10), ha='center')
        ax.set_xlim(200, 0)
        plt.show()
        plt.close('all')
    # Integrate Peaks - for each compound, integrate the simulated signal of
    # just the peaks belonging to that compound
    est_signals = {}
    real_signals = {}
    for cpd in all_cpds:
        ran = np.full((data.shape[0],), False)
        include = np.isin(shifts, assignments[cpd])
        idxs = assignment_map[cpd]
        print([ixc(sh) for sh in idxs])
        for sh in idxs:
            ran = ran | (ppm <= ixc(sh) + 0.45) & (ppm >= ixc(sh) - 0.45)
        cpdsignal = copy.deepcopy(data)
        cpdsignal[~ran] = 0.
        # sim = ng.analysis.linesh.sim_NDregion( # Simulate spectrum
        #     data.shape, ['v'], params[include], amps[include])
        # areacurve = ng.analysis.integration.integrate(sim, uc, (0, 200))
        areashape = ng.analysis.integration.integrate(cpdsignal, uc, (0, 200))
        # est_signals[cpd] = areacurve[-1]
        real_signals[cpd] = areashape[-1]
    return None, real_signals # est_signals,
