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
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors
import openpyxl as xl

from curveshapes import logi_flex as fx
from curveshapes import ddx_logi_flex as dfdx
from synchronize import synchronizers
from trajectories import call_fit, fit_trajectories

SCDIR = os.path.dirname(__file__)   # location of script
BOLD = xl.styles.Font(bold=True)

tracked_reactions = [
    "ID_391", # Enolase
    "ID_53",  # PFOR
    "ID_280", # Ac Kinase
    "ID_326", # WLP
    "RNF-Complex", # RNF
    "ID_336", # ALT
    "RXN-19534",  # Ox Leu
    "ICCoA-DHG-EB",  # Red Leu
    "ID_314",  # Pro Red
    "ID_383",  # EtOH dehydrogenase
    "BUK",     # Butyrate kinase
    "Sec_h2",   # h2 evolution
    "ATP_sink", # ATP objective
    "ATPsynth4_1", #ATP synthase
    "ID_575", #GDH
    # "Ex_biomass",
]

tracked_metabolites = [
    "alaL", #L-alanine
    "gluL", #L-glutamate
    "atp",  #ATP
    "pyr",  #Pyruvate
    "nh3",  #Ammonia
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
}

def find_minmax_bounds(x, bounds):
    """Estimate min/max bounds of function using coefficient error.
    Uses an asinine and wasteful methodology.
    """
    all_bounds = []
    all_signal = []
    options = zip(*bounds)
    combinations = list(itertools.product(*options))
    for coeffset in combinations:
        all_bounds.append(dfdx(x, *coeffset))
        all_signal.append(fx(x, *coeffset))
    return min(all_bounds), max(all_bounds), min(all_signal), max(all_signal)

def set_bounds(model, rid, lower=0., upper=1000., update=True):
    if update:
        rxn = model.reactions.get_by_id(rid)
        rxn.bounds = (lower, upper)

def update_uptake_bounds(model, t, met, coeffs, coeff_bounds, update=True):
    """Update the uptake rate of exchange reaction at timepoint t."""

    # Calculate logistic solutions with flexible number of parameters
    signal = fx(t, *coeffs)
    exch = dfdx(t, *coeffs)
    exch_l, exch_u, signal_l, signal_u = find_minmax_bounds(t, coeff_bounds + [coeffs])
    # Reverse direction for secretion reactions
    if met in ['glc', 'proL', 'leuL', 'valL']:
        rid = 'Ex_' + met
        exch *= -1
        tmp = exch_l
        exch_l = -1*exch_u
        exch_u = -1*tmp
        if met in ['leuL', 'valL']:
            set_bounds(model, rid, update=update)
        else:
            set_bounds(model, rid, lower=exch_l, upper=exch_u, update=update) # lower=exch_l, upper=exch_u
    elif met in ['ival', 'isocap', 'isobuta', '5av', 'but']:
        rid = 'Sec_' + met
        set_bounds(model, rid, lower=exch, upper=exch, update=update)
    else:
        rid = 'Sec_' + met
        set_bounds(model, rid, lower=exch_l, update=update) # lower=exch_l,
    return exch, exch_l, exch_u, signal, signal_l, signal_u

def areaplot(df, dl, du, t_max=48, ylabel='flux (mol/gDW/h)'):
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

def dfba_main(all_params, traced, fba_method=lsol, prop_ids=["alaL", "gluL"],
              modelfile=f"{SCDIR}/../../data/icdf843.json", t_max=48,
              resolution=1, obj="ATP_sink", dry_run=False):
    """Main dFBA function. Computes successive FBA solutions and plots the concentration, uptake rates, and tracked reaction fluxes.

       Arguments:
       specsheet -- path to tsv sheet containing logistic parameters
       traced -- list of reaction ids to be traced over the timecourse

       Keyword arguments:
       modelfile -- location of metabolic model
       t_max -- end timepoint in hours (default 48)
       resolution -- number of hours per dFBA sample (default 1)
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

    # pool_rxn = cb.Reaction(
    #     'NADH_pool', name='Electron pooling', lower_bound=0, upper_bound=0
    # )
    # model.add_reaction(pool_rxn)
    # pool_rxn.add_metabolites({'nadh_c':-1, 'nad_c':1, 'h_c':1})

    init_cnc = dict()
    for rxn in model.reactions:
        if rxn.id.startswith('Ex_') and rxn.id.endswith('L') \
                or rxn.id in ['Ex_gly', 'Ex_his'] \
                and rxn.id != 'Ex_valL':
            init_cnc[rxn.id] = rxn.upper_bound
            rxn.upper_bound *= 0.03
    model.reactions.Ex_glc.upper_bound=0

    params = {met: all_params[0][name] for name, met in metmap.items()}
    par_lb = {met: all_params[1][name] for name, met in metmap.items()}
    par_ub = {met: all_params[2][name] for name, met in metmap.items()}

    timecourse = list(np.linspace(0, t_max, num=t_max*resolution+1))
    mets = [k for k in params]

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
        np.zeros(shape=(len(timecourse), len(params))),
        index=timecourse,
        columns=[model.metabolites.get_by_id(k + '_c').name for k in mets]
    )
    results = [datashape]
    for i in range(5):
        results.append(datashape.copy(deep=True))

    rxnflux = pd.DataFrame(
        np.zeros(shape=(len(timecourse), len(traced))),
        index=timecourse,
        columns=[model.reactions.get_by_id(rid).name for rid in traced]
    )
    if fva_run:
        rxnf_ub = rxnflux.copy(deep=True)
        rxnf_lb = rxnflux.copy(deep=True)

    print('dFBA log: simulation output begin.')

    for i, t in enumerate(timecourse):
        if t < params['isocap'][2]:
            model.reactions.ID_326.upper_bound=0
        else:
            model.reactions.ID_326.upper_bound=1000

        for j, met in enumerate(mets):
            coeffs = params[met]
            coeff_bounds = [par_lb[met], par_ub[met]]
            fbounds = update_uptake_bounds(model, t, met, coeffs, coeff_bounds)
            mk = model.metabolites.get_by_id(met + '_c').name
            for k, df in enumerate(results):
                df.at[t, mk] = fbounds[k]
            # if met == 'ival':
            #     t3 = update_uptake_bounds(model, t-3, met, coeffs, coeff_bounds,
            #                               update=False)
            #     set_bounds(model, "NADH_pool", lower=1.72*(fbounds[1]-t3[2]), upper=1.72*(fbounds[2]-t3[1]))

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
        for rid in traced:
            rname = model.reactions.get_by_id(rid).name
            rxnflux.at[t, rname] = sol.fluxes[rid]
            if fva_run:
                rxnf_ub.at[t, rname] = sol_v.at[rid, 'maximum']
                rxnf_lb.at[t, rname] = sol_v.at[rid, 'minimum']

        if i in [0, 6, 8, 12, 24, 36, 48]:
            print("time = " + str(i))
            print(model.summary(solution=sol))
            for met in propmets:
                print(met.summary(solution=sol))

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
        if fva_run:
            print(f"t = {t} simulation complete")
    # Adjust flux solution directionality in output for visualization
    reverse_ids = ["ID_383", "ID_336", "ID_391", "HydEB"]
    reverse_nms = [model.reactions.get_by_id(rid).name for rid in reverse_ids]

    for rname in reverse_nms:
        if rname in rxnflux:
            rxnflux[rname] *= -1
            if fva_run:
                rxnf_ub[rname] *= -1
                rxnf_lb[rname] *= -1

    print(f'Complete after {i} cycles ({t} hours). Final flux: {sol.fluxes["Ex_biomass"]}')
    areaplot(results[0], results[1], results[2], ylabel='flux (mM/h)')
    plt.show()
    areaplot(results[3], results[4], results[5], ylabel='normalized signal')
    plt.show()
    niceplot(rxnflux)
    plt.show()
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

methods = {
    'fba': lambda m: m.optimize(),
    'fva': fva,
    'loopless': lsol,
    'pfba': pfba
}

tracked = tracked_reactions

dfba_kwargs = {
    'modelfile': f'{SCDIR}/../../data/icdf843.json',
    't_max': 36,
    'resolution': 1,
    'obj': 'ATP_sink',
    # 'obj': 'Ex_biomass',
}

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 1:
        print("Please provide a FBA method to use. Options: " +
              ", ".join(methods))
    else:
        method = args.pop(0)
        if method not in methods:
            print(f"Unsupported FBA method {method}. Please select one of "
                  + ", ".join(methods))
        elif len(args) == 0:
            dfba_main(call_fit(), tracked, fba_method=methods[method], **dfba_kwargs)
        else:
            dirs = [fp.rstrip('/') for fp in args]
            params = [{}, {}, {}]
            for i, tscale in enumerate(synchronizers(dirs, stretch=True, plot=True)):
                filepath = dirs[i]
                basename = os.path.basename(filepath)
                curves = fit_trajectories(filepath, tscale=tscale, plot=True)
                for j, ele in enumerate(curves):
                    params[j].update(ele)
            dfba_main(
                params, tracked, fba_method=methods[method], **dfba_kwargs,
                prop_ids=tracked_metabolites
            )
