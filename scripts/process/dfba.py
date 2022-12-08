#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  22 16:31:12 2022

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Run dynamic flux balance analysis (dFBA) with NMR constraints
 - Use python 3.8+ for best results
 - See /nmr-cdiff/venv/requirements.txt for dependencies

Copyright 2022 Massachusetts Host-Microbiome Center

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
import itertools
import os
import sys

import cobra as cb
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.flux_analysis.loopless import loopless_solution as lsol
from cobra.flux_analysis import pfba
from cobra.flux_analysis.variability import find_blocked_reactions as fbr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors
import openpyxl as xl

from curveshapes import logi_flex as fx
from curveshapes import ddx_logi_flex as dfdx
from curveshapes import err_logi_flex as del_f
from curveshapes import err_ddx_logi_flex as del_df
from synchronize import synchronizers
from trajectories import call_fit, fit_trajectories

SCDIR = os.path.dirname(__file__)   # location of script
BOLD = xl.styles.Font(bold=True)

tracked_reactions = [
    "ID_233", # 1. Phosphoglycerate kinase
    "ID_53",  # 2. PFOR
    "ID_280", # 3. Ac Kinase
    "ID_326", # 4. WLP
    "RNF-Complex", # 5. RNF
    "ID_336", # 6. ALT
    "ID_366",  # 7. isovalerate kinase
    "ICCoA-DHG-EB",  # 8. Red Leu
    "ID_314",  # 9. Pro Red
    "ID_383",  # 10. EtOH dehydrogenase
    "BUK",     # 11. Butyrate kinase
    "HydEB",   # 12. electron-bifurcating hydrogenase
    "ATP_sink", # 13. ATP objective
    "ATPsynth4_1", # 14. ATP synthase
    "ID_575", # 15. GDH
    "ID_469", # 16. Cys gamma-lyase
    "ID_146", # 17. Ox Ile
    "ID_321", # 18. Ox Val
    "ID_252", # 19. pyruvate kinase
    "ID_407", # UDP
    "ID_582", # CDP
    "ID_1",   # IDP
    "Ex_biomass",
]

tracked_metabolites = [
    "alaL", #L-alanine
    "gluL", #L-glutamate
    "atp",  #ATP
    "pyr",  #Pyruvate
    "nh3",  #Ammonia
    "feroxred", #Ferredoxin
]

metmap = {
    'Glucose': 'glc',
    'Acetate': 'ac',
    'Alanine': 'alaL',
    'Ethanol': 'eto',
    'Butyrate': 'but',
    'Proline': 'proL',
    'Leucine': 'leuL',
    'Isovalerate': 'ival',
    'Isocaproate': 'isocap',
    'Valine': 'valL',
    'Isobutyrate': 'isobuta',
    'Isoleucine': 'ileL',
    '2-methylbutyrate': '2mbut',
    'Acetate13C': 'acoa',
    'AcetateNA': 'gly',
}

stretch = False

def find_minmax_bounds(x, popts, perrs):
    """Estimate min/max bounds of function using coefficient error.
    Arguments:
    x -- time series vector (array) or an individual timepoint (int or float)
    popts -- optimal logistic coefficients
    perrs -- logistic coefficient SEMs

    Returns 95% confidence interval of exchange fluxes and estimated
    concentrations for each element in the time series vector.
    """
    f_vals = fx(x, *popts)
    df_vals = dfdx(x, *popts)
    f_errs = del_f(x, *popts, *perrs)
    df_errs = del_df(x, *popts, *perrs)
    exch_l = df_vals - 1.96*df_errs
    exch_u = df_vals + 1.96*df_errs
    signal_l = f_vals - 1.96*f_errs
    signal_u = f_vals + 1.96*f_errs
    return exch_l, exch_u, signal_l, signal_u

def set_bounds(model, rid, lower=0., upper=1000., update=True):
    """Set upper and lower flux bounds for a reaction.

    Parameters:
    model -- COBRA model
    rid -- the ID of the reaction to bound
    lower -- lower bound value to set (default: 0)
    upper -- upper bound value to set (default: 1000)
    update -- whether to set the bounds (default: True)
    """
    if update:
        rxn = model.reactions.get_by_id(rid)
        rxn.bounds = (lower, upper)

def update_uptake_bounds(model, t, met, popts, perrs, update=True):
    """Calculate and set exchange reaction bounds at timepoint t.

    Parameters:
    model -- COBRA model
    t -- the individual timepoint (int or float)
    met -- the metabolite for which the exchange bounds are to be changed
    popts -- optimal logistic coefficients
    perrs -- logistic coefficient SEMs
    update -- whether to set the bounds (default: True)

    Returns optimal values and 95% confidence interval of exchange fluxes and
    estimated concentrations at the timepoint. If update==True, also sets the
    exchange reaction bounds to the 95% confidence interval limits.
    """
    # Calculate logistic solutions with flexible number of parameters
    signal = fx(t, *popts)
    exch = dfdx(t, *popts)
    exch_l, exch_u, signal_l, signal_u = find_minmax_bounds(t, popts, perrs)
    # Reverse direction for secretion reactions
    if met in ['glc', 'proL', 'leuL', 'valL', 'ileL']:
        rid = 'Ex_' + met
        exch *= -1
        tmp = exch_l
        exch_l = -1*exch_u
        exch_u = -1*tmp
        if met in ['leuL']:
            set_bounds(model, rid, update=update)
        else:
            set_bounds(model, rid, lower=exch_l, upper=exch_u, update=update) # lower=exch_l, upper=exch_u
    elif met in ['but']:
        set_bounds(model, 'Sec_but', lower=exch_l, upper=exch_u, update=update)
        set_bounds(model, 'ID_326', lower=exch_l, upper=exch_u, update=update)
    elif met in ['gly', 'acoa']:
        pass
    else:
        rid = 'Sec_' + met
        lb = max(exch_l, 0)
        ub = max(exch_u, 0)
        set_bounds(model, rid, lower=lb, upper=ub, update=update) # lower=exch_l,
    return exch, exch_l, exch_u, signal, signal_l, signal_u

def areaplot(df, dl, du, t_max=48, ylabel='flux (mol/gDW/h)'):
    """Plot lines with shaded confidence interval.

    Parameters:
    df -- dataframe of vectors corresponding to optimal values
    dl -- dataframe of vectors corresponding to lower bounds
    du -- dataframe of vectors corresponding to upper bounds
    t_max -- maximum value of x-axis
    ylabel -- y-axis label
    """
    f = plt.figure(figsize=(14, 10), )#fontsize=30,
    ax = plt.axes(
        xticks=(range(0, t_max+1, 12)),
        xlim=(0, t_max),
        ylim=(0, 1.05*max(ele.to_numpy().max() for ele in [df, dl, du])),
    )
    for met, ser in df.iteritems():
        ser.index.name = 'index'
        p = ax.plot(ser.index.to_numpy(), ser.to_numpy(), '-', label=met, lw=5)
        color = matplotlib.colors.to_rgb(p[0].get_color())
        color = color + (0.2,)
        ax.fill_between(dl.index.to_numpy(), dl[met].to_numpy(), du[met].to_numpy(), color=color)
    ax.set_xlabel('time (h)', fontsize=30, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=30, fontweight='bold')
    plt.xticks(ax.get_xticks(), weight='bold')
    plt.yticks(ax.get_yticks(), weight='bold')
    plt.show()

def niceplot(df, t_max=48, ylabel='flux (mol/gDW/h)'):
    """Plot lines.

    Parameters:
    df -- dataframe of data to plot
    t_max -- maximum value of x-axis
    ylabel -- y-axis label
    """
    ax = df.plot(
        figsize=(14, 10),
        xticks=(range(0, t_max+1, 12)),
        xlim=(0, t_max),
        ylim=(0, None),
        fontsize=30,
        lw=5,
    )
    ax.set_xlabel('time (h)', fontsize=30, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=30, fontweight='bold')
    plt.xticks(ax.get_xticks(), weight='bold')
    plt.yticks(ax.get_yticks(), weight='bold')
    plt.title("dFBA (10/h), ATP, WLP > 8h")
    plt.show()

def dfba_main(params_opt, params_err, tracked_rxns, fba_method=lsol, prop_ids=["alaL", "gluL"],
              modelfile=f"{SCDIR}/../../data/icdf843.json", t_max=48,
              resolution=1, obj="ATP_sink", dry_run=False):
    """Main function to compute dFBA solutions.
    Computes successive static FBA solutions and plots the estimated metabolite
    concentrations, uptake rates, and tracked reaction fluxes.

    Parameters:
    params_opt -- dictionary mapping metabolite names to optimal logistic
                  coefficients
    params_err -- dictionary mapping metabolite names to standard errors of
                  logistic coefficients
    tracked_rxns -- reactions to tracked and written to fluxes.xlsx
    fba_method -- FBA method to use for static solutions (function, default:
                  cobra.flux_analysis.loopless.loopless_solution)
    modelfile -- location of metabolic model
    t_max -- end timepoint in hours (default 48)
    resolution -- number of static solutions per hour (default 1)
    obj -- reaction ID of objective function (default ATP_sink)
    dry_run -- True to avoid writing output to file (default False)
    """

    print(f"""dFBA log: Begin dFBA analysis with sheet, endpoint
          {t_max} hours, and resolution {resolution}.""")
    print(f"Using method {fba_method}.")
    if dry_run:
        print("Dry run, will not write results.")
    # Load metabolic model and logistic fit specs #
    print('dFBA log: loading model and specsheet...')
    if fba_method == fva:
        fva_run = True
    else:
        fva_run = False

    model = cb.io.load_json_model(modelfile)
    model.objective = obj
    model.reactions.get_by_id(obj).upper_bound=1000

    init_cnc = dict()
    for rxn in model.reactions:
        if ( rxn.id.startswith('Ex_') and rxn.id.endswith('L') \
                or rxn.id in ['Ex_gly', 'Ex_his'] ):
            init_cnc[rxn.id] = rxn.upper_bound
            rxn.upper_bound *= 0.03
        if rxn.id in ['Ex_valL', 'Ex_ileL']:
            rxn.upper_bound=0
    model.reactions.Ex_glc.upper_bound=0
    model.reactions.Ex_thrL.upper_bound = 0
    model.reactions.Ex_cysL.upper_bound = 1000

    popts = {m: params_opt[n] for n, m in metmap.items() if n in params_opt}
    perrs = {m: params_err[n] for n, m in metmap.items() if n in params_err}

    for met in ("valL", "ileL", "isobuta", "2mbut"):
        for dic in popts, perrs:
            if met in ("valL", "isobuta"):
                dic[met][2] = dic["isocap"][2]
            else:
                dic[met][2] = dic["proL"][2]

    timecourse = list(np.linspace(0, t_max, num=t_max*resolution+1))
    mets = list(popts)

    propdata = dict()
    propmets = []
    for met_id in prop_ids:
        if model.metabolites.has_id(met_id + "_c"):
            propdata[met_id] = {
                'rxns_in': set([]),
                'rxns_out': set([]),
                'data_in': [],
                'data_out': [],
            }
            propmets.append(model.metabolites.get_by_id(met_id + "_c"))
            if fva_run:
                propdata[met_id]['data_in_lb'] = []
                propdata[met_id]['data_in_ub'] = []
                propdata[met_id]['data_out_lb'] = []
                propdata[met_id]['data_out_ub'] = []
        else:
            print(f"'{met_id}_c' is not a valid metabolite.")

    # Initialize results dataframes
    datashape = pd.DataFrame(
        np.zeros(shape=(len(timecourse), len(popts))),
        index=timecourse,
        columns=[model.metabolites.get_by_id(k + '_c').name for k in mets]
    )
    results = [datashape]
    for i in range(5):
        results.append(datashape.copy(deep=True))

    rxnflux = pd.DataFrame(
        np.zeros(shape=(len(timecourse), len(tracked_rxns))),
        index=timecourse,
        columns=[model.reactions.get_by_id(rid).name for rid in tracked_rxns]
    )
    if fva_run:
        rxnf_ub = rxnflux.copy(deep=True)
        rxnf_lb = rxnflux.copy(deep=True)

    print('dFBA log: simulation output begin.')

    for i, t in enumerate(timecourse):
        for j, met in enumerate(mets):
            ## Calculate and set exchange constraints ##
            fbounds = update_uptake_bounds(model, t, met, popts[met], perrs[met])
            mk = model.metabolites.get_by_id(met + '_c').name

            ## Record exchange constraints ##
            # if met != 'glc':
            for k, df in enumerate(results):
                df.at[t, mk] = fbounds[k]

        #     ## Calculate glucose uptake constraint ##
        #     if met in ['ac', 'eto', 'but', 'alaL']:
        #         for k, df in enumerate(results):
        #             if met == 'but':
        #                 df.at[t, 'beta-D-glucose'] += fbounds[k]
        #             else:
        #                 df.at[t, 'beta-D-glucose'] += 0.5*fbounds[k]
        #
        # ## Set glucose uptake constraint ##
        # set_bounds(model, "Ex_glc", lower=results[1].at[t, 'beta-D-glucose'],
        #            upper=results[2].at[t, 'beta-D-glucose'], update=True)

        ## Add natural abundance acetate allowance to secretion constraint ##
        lower, upper = model.reactions.Sec_ac.bounds
        set_bounds(model, "Sec_ac", lower=lower+results[1].at[t, 'glycine'],
                   upper=upper+results[2].at[t, 'glycine'], update=True)

        ## Get flux solution(s) and populate arrays ##
        if fva_run:
            sol = lsol(model)
            sol_v = fba_method(model, fraction_of_optimum=0.995, loopless=True)
        else:
            sol = fba_method(model)
        if i % 10 == 0:
            print(f'dFBA log: Time = {t}  (cycle {i+1}) \tFBA solution: {sol.fluxes[obj]}')
        if sol.fluxes[obj] < 0.0001:
            print(f'dFBA log: infeasible solution on cycle {i}.')
        for rid in tracked_rxns:
            rname = model.reactions.get_by_id(rid).name
            rxnflux.at[t, rname] = sol.fluxes[rid]
            if fva_run:
                rxnf_ub.at[t, rname] = sol_v.at[rid, 'maximum']
                rxnf_lb.at[t, rname] = sol_v.at[rid, 'minimum']

        ## Display incremental solutions ##
        if t in [0, 6, 8, 12, 24, 36, 48]:
            print("time = " + str(i))
            print(model.summary(solution=sol))
            for met in propmets:
                print(met.summary(solution=sol))

        ## Record flux contributions for tracked metabolites ##
        for ac in propmets:
            met_id = ac.id[:-2]
            fldic = {}
            for dr in ('in', 'out'):
                labels = [f"data_{dr}"]
                if fva_run:
                    labels.extend([f"data_{dr}_lb", f"data_{dr}_ub"])
                fldic.update({k: dict() for k in labels})
                for rxn in ac.reactions:
                    if dr == 'out':
                        side = (rxn.reactants, rxn.products)
                    else:
                        side = (rxn.products, rxn.reactants)
                    if (ac in side[0] and sol.fluxes[rxn.id] > 1E-5) \
                       or (ac in side[1] and sol.fluxes[rxn.id] < -1E-5):
                        fldic[f"data_{dr}"][rxn.id] = rxn.metabolites[ac] * sol.fluxes[rxn.id]
                        if fva_run:
                            fldic[f"data_{dr}_lb"][rxn.id] = rxn.metabolites[ac] * sol_v.at[rxn.id, 'minimum']
                            fldic[f"data_{dr}_ub"][rxn.id] = rxn.metabolites[ac] * sol_v.at[rxn.id, 'maximum']
                propdata[met_id][f'rxns_{dr}'].update(fldic[f"data_{dr}"].keys())
                if fva_run:
                    propdata[met_id][f'rxns_{dr}'].update(fldic[f"data_{dr}_lb"].keys())
                    propdata[met_id][f'rxns_{dr}'].update(fldic[f"data_{dr}_ub"].keys())
                for lab in labels:
                    propdata[met_id][lab].append(fldic[lab])

        # FVA takes longer, print at end of each simulation
        if fva_run:
            print(f"t = {t} simulation complete")

    ## Adjust flux solution directionality in output for visualization ##
    reverse_ids = ["ID_383", "ID_336", "ID_391", "HydEB"]
    reverse_nms = [model.reactions.get_by_id(rid).name for rid in reverse_ids]
    for rname in reverse_nms:
        if rname in rxnflux:
            rxnflux[rname] *= -1
            if fva_run:
                rxnf_ub[rname] *= -1
                rxnf_lb[rname] *= -1

    ## Plot run results ##
    print(f'Complete after {i} cycles ({t} hours). Final flux: {sol.fluxes["Ex_biomass"]}')
    areaplot(results[0], results[1], results[2], ylabel='flux (mM/h)')
    plt.show()
    areaplot(results[3], results[4], results[5], ylabel='normalized signal')
    plt.show()
    niceplot(rxnflux)
    plt.show()

    ## Write run results ##
    if not dry_run:
        writer = pd.ExcelWriter(f'{SCDIR}/../../data/fluxes.xlsx', engine='openpyxl')
        results[0].to_excel(writer, sheet_name='data', engine='openpyxl')
        results[1].to_excel(writer, sheet_name='lbdev', engine='openpyxl')
        results[2].to_excel(writer, sheet_name='ubdev', engine='openpyxl')
        results[3].to_excel(writer, sheet_name='signal', engine='openpyxl')
        results[4].to_excel(writer, sheet_name='slbdev', engine='openpyxl')
        results[5].to_excel(writer, sheet_name='subdev', engine='openpyxl')
        rxnflux.to_excel(writer, sheet_name='fluxes', engine='openpyxl')
        if fva_run:
            rxnf_ub.to_excel(writer, sheet_name='fluxub', engine='openpyxl')
            rxnf_lb.to_excel(writer, sheet_name='fluxlb', engine='openpyxl')
        writer.close()

    if not dry_run:
        writer = pd.ExcelWriter(f'{SCDIR}/../../data/met_fluxes.xlsx', engine='openpyxl')
        wb=writer.book
        ws=wb.create_sheet("MetFluxes")
        writer.sheets["MetFluxes"] = ws
        scol = 0

    for ac in propmets:
        met_id = ac.id[:-2]
        print(f"Flux slices for metabolite {ac.name} ({ac.id}).")
        def make_df(cols):
            return pd.DataFrame(data=np.zeros((len(timecourse), len(cols))),
                                index=timecourse, columns=cols)
        fldfs = {}
        rnms = {}
        rids = {}
        names = ["Fluxes", "Flux Upper Bound", "Fluxes Lower Bound"]

        for dr in ('in', 'out'):
            rids[dr] = list(propdata[met_id][f'rxns_{dr}'])
            rnms[dr] = [model.reactions.get_by_id(rid).name for rid in rids[dr]]
            labels = [f'data_{dr}']
            if fva_run:
                labels.extend([f'data_{dr}_lb', f'data_{dr}_ub'])
            fldfs.update({k: make_df(rids[dr]) for k in labels})
            for i, t in enumerate(timecourse):
                for rxn in rids[dr]:
                    for k in labels:
                        fldfs[k].at[t, rxn] = propdata[met_id][k][i].get(rxn, 0)
            for i, lab in enumerate(labels):
                print(names[i])
                print(f"\t{dr}flux percentages:")
                v = fldfs[lab]
                totalflux = v.to_numpy().sum()
                cuts = v.sum(axis=0)
                for i, rxn in enumerate(rids[dr]):
                    print("\t", rxn, rnms[dr][i], f"{cuts.at[rxn]/totalflux*100:.2f}")
                print(f"\tAbsolute {dr}flux:")
                for i, rxn in enumerate(rids[dr]):
                    print("\t", rxn, rnms[dr][i], f"{cuts.at[rxn]}")
            if not dry_run:
                df = fldfs[f'data_{dr}']
                ecol = scol + df.shape[1]
                if scol == 0:
                    w_index = True
                    ecol += 1
                    ws.cell(row=1, column=2, value="Time")
                else:
                    w_index = False
                df.to_excel(writer, sheet_name="MetFluxes", startrow=2,
                            startcol=scol, index=w_index)
                c = ws.cell(row=1, column=scol+1, value=f"{ac.name} {dr}flux")
                c.font = BOLD
                ws.merge_cells(start_row=1, start_column=scol+1, end_row=1,
                               end_column=ecol)
                colnames = list(df.columns)
                for i, j in enumerate(range(ecol-df.shape[1], ecol)):
                    rxname = model.reactions.get_by_id(colnames[i]).name
                    c = ws.cell(row=2, column=j+1, value=rxname)
                    c.font = BOLD
                scol = ecol
                with open(f'{SCDIR}/../../data/{met_id}flux_{dr}.txt', 'w') as wf:
                    wf.write('time\t' + '\t'.join(rids[dr]) + '\n')
                    wf.write('time\t' + '\t'.join(rnms[dr]) + '\n')
                    for i, t in enumerate(timecourse):
                        wf.write(str(t))
                        for rxn in rids[dr]:
                            val = df.at[t, rxn]
                            wf.write('\t' + str(val))
                        wf.write('\n')
    if not dry_run:
        writer.close()

METHODS = {
    'fba': lambda m: m.optimize(),
    'fva': fva,
    'loopless': lsol,
    'pfba': pfba
}

DFBA_KWARGS = {
    'modelfile': f'{SCDIR}/../../data/icdf843.json',
    't_max': 36,
    'resolution': 5,
    'obj': 'ATP_sink',
    # 'obj': 'Ex_biomass',
}

if __name__ == "__main__":
    """Handle command line call"""
    args = sys.argv[1:]
    if len(args) < 1:
        print("Please provide a FBA method to use. Options: " +
              ", ".join(METHODS))
    else:
        method = args.pop(0)
        if method not in METHODS:
            print(f"Unsupported FBA method {method}. Please select one of "
                  + ", ".join(METHODS))
        elif len(args) == 0:
            dfba_main(call_fit(), tracked_reactions, fba_method=METHODS[method],
                      **DFBA_KWARGS)
        else:
            dirs = [fp.rstrip('/') for fp in args]
            params = [{}, {}]
            for i, tscale in enumerate(synchronizers(dirs, stretch=stretch, \
                    plot=False)):
                filepath = dirs[i]
                basename = os.path.basename(filepath)
                curves = fit_trajectories(filepath, tscale=tscale, plot=True)
                for j, ele in enumerate(curves):
                    params[j].update(ele)
            dfba_main(
                *params, tracked_reactions, fba_method=METHODS[method],
                **DFBA_KWARGS, prop_ids=tracked_metabolites
            )
