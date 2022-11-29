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
from collections import defaultdict
import os
import sys

import cobra as cb
from matplotlib import pyplot as plt
import matplotlib.colors
import numpy as np
import openpyxl as xl
import pandas as pd

from curveshapes import LogisticSet
from synchronize import synchronizers
from trajectories import fit_trajectories
from get_color import get_cmap

SCDIR = os.path.dirname(__file__)   # location of script
BOLD = xl.styles.Font(bold=True)
MODEL_PATH = f'{SCDIR}/../../data/icdf843.json'
METHODS = {
    'fba': lambda m: m.optimize(),
    'fva': cb.flux_analysis.flux_variability_analysis,
    'loopless': cb.flux_analysis.loopless.loopless_solution,
    'pfba': cb.flux_analysis.pfba,
}

## Objective function for dFBA simulations
objective_function = 'ATP_sink'
# objective_function = 'Ex_biomass'
# objective_function = 'Sec_exopoly'

## Arguments to dFBA function. Edit these to modify behavior.
dFBA_kwargs = {
    'modelfile': MODEL_PATH, # string filepath to model
    't_max': 36, # end of simulated timecourse, in hours
    'resolution': 5, # number of solutions per hour in simulated timecourse
    'obj': objective_function, # Objective function, change above
}

## Identifiers of reactions to track
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
    "ID_252", # 19. pyruvate kinase (ADP)
    "ID_407", # pyruvate kinase (UDP)
    "ID_582", # pyruvate kinase (CDP)
    "ID_1",   # pyruvate kinase (IDP)
    "Ex_biomass",
    "Sec_exopoly",
]

## Identifiers of metabolites to track (i.e. compute reaction flux fractions)
tracked_metabolites = [
    "alaL", #L-alanine
    "gluL", #L-glutamate
    "atp",  #ATP
    "pyr",  #Pyruvate
    "nh3",  #Ammonia
    "datp",
    "fru16bp",
]

## Map names to metabolite identifiers, add to this as new datasets are included
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

namemap = {
    'Glucose': 'beta-D-glucose',
    'Acetate': 'acetate',
    'Alanine': 'L-alanine',
    'Ethanol': 'ethanol',
    'Butyrate': 'Butyrate',
    'Proline': 'L-proline',
    'Leucine': 'L-leucine',
    'Isovalerate': 'Isovalerate',
    'Isocaproate': 'Isocaproate',
    'Valine': 'L-valine',
    'Isobutyrate': 'isobutyric acid',
    'Isoleucine': 'L-isoleucine',
    '2-methylbutyrate': '2-Methylbutyric acid',
    'Acetate13C': 'acetyl-CoA',
    'AcetateNA': 'glycine',
}

## Whether to stretch (True) or just shift (False) timecourse during normalization
stretch = False

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

def reverse_flux(flux, lb, ub):
    """Reverse flux direction and swap upper/lower bounds."""
    return -1*flux, -1*ub, -1*lb

def update_uptake_bounds(model, t, met, params, update=True):
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
    signal, serr, exch, exerr = params.get_sol(t)
    exch_l, exch_u, signal_l, signal_u = params.get_bounds(t)
    # Reverse direction and set bounds for secretion reactions
    if met in ['glc', 'proL', 'leuL', 'valL', 'ileL']:
        rid = 'Ex_' + met
        exch, exch_l, exch_u = reverse_flux(exch, exch_l, exch_u)
        if met in ['leuL']:
            set_bounds(model, rid, update=update) # leave leucine unbounded
        else:
            set_bounds(model, rid, lower=exch_l, upper=exch_u, update=update)
    # Update Wood-Ljungdahl Pathway bounds with butyrate
    elif met in ['but']:
        set_bounds(model, 'Sec_but', lower=exch_l, upper=exch_u, update=update)
        set_bounds(model, 'ID_326', lower=exch_l, upper=exch_u, update=update)
    # Ignore data from 1H spectra
    elif met in ['gly', 'acoa']:
        pass
    # Update bounds for products
    else:
        rid = 'Sec_' + met
        lb = max(exch_l, 0) # do not allow reverse flux
        ub = max(exch_u, 0)
        set_bounds(model, rid, lower=lb, upper=ub, update=update)
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
        ax.fill_between(
            dl.index.to_numpy(), 
            dl[met].to_numpy(), 
            du[met].to_numpy(), 
            color=color
        )
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

def areaplot2(t, params):
    """Plot lines with shaded confidence interval.

    Parameters:
    df -- dataframe of vectors corresponding to optimal values
    dl -- dataframe of vectors corresponding to lower bounds
    du -- dataframe of vectors corresponding to upper bounds
    t_max -- maximum value of x-axis
    ylabel -- y-axis label
    """
    substrate_map = {
        'Glucose': [
            'Glucose',
            'Acetate',
            'Alanine',
            'Ethanol',
            'Butyrate',
        ],
        'Leucine': [
            'Leucine',
            'Isovalerate',
            'Isocaproate',
        ],
        'Proline': [
            'Proline',
            '5-aminovalerate',
        ],
    }
    cmap = get_cmap()
    for substrate, products in substrate_map.items():
        f = plt.figure(figsize=(2.5, 2), constrained_layout=True)
        ax = plt.axes()
        for met in products:
            ser = params[met].get_sol(t)[0]
            _, _, lb, ub = params[met].get_bounds(t)
            color = cmap[met]
            ax.plot(t, ser, '-', label=met, lw=2, c=color)
            color = color + (0.2,)
            ax.fill_between(t, lb, ub, color=color)
        afont = {'fontname': 'Arial', 'size': 7}
        ax.set_xlabel("Time (h)", **afont)
        ax.set_ylabel("Estimated Concentration (mM)", **afont)
        ax.set_xlim((0, 36))
        if substrate == "Glucose":
            ax.set_yticks([0, 10, 20, 30])
            ax.set_yticks(list(range(-2, 31)), minor=True)
        elif substrate == "Proline":
            ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
            ax.set_yticks(list(np.arange(-0.5, 8, 0.5)), minor=True)
        elif substrate == "Leucine":
            ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
            ax.set_yticks(list(np.arange(-0.5, 9, 0.5)), minor=True)
        ax.set_xticks([0, 12, 24, 36])
        ax.set_xticks(list(range(37)), minor=True)
        ax.set_xticklabels(ax.get_xticks(), **afont)
        ax.set_yticklabels(ax.get_yticks(), **afont)
        ax.xaxis.set_tick_params(width=0.5)
        ax.yaxis.set_tick_params(width=0.5)
        plt.setp(ax.spines.values(), linewidth=0.5)
        plt.show()

def dfba_main(params, tracked_rxns, fba_method=METHODS['fba'], fva_run=False, 
              tracked_mets=["alaL", "gluL"], modelfile=MODEL_PATH,
              t_max=48, resolution=1, obj="ATP_sink", dry_run=False):
    """Main function to compute dFBA solutions.
    Computes successive static FBA solutions and plots the estimated metabolite
    concentrations, uptake rates, and tracked reaction fluxes.

    Parameters:
    params -- dictionary mapping metabolite names to LogisticSet objects
            containing optimal logistic coefficients and errors
    tracked_rxns -- reactions to tracked and written to fluxes.xlsx
    fba_method -- FBA method to use for static solutions (function, default:
            cobra.flux_analysis.loopless.loopless_solution)
    fva_run -- whether to also compute an FVA solution
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
    model = cb.io.load_json_model(modelfile)
    model.objective = obj
    model.reactions.get_by_id(obj).upper_bound=1000

    # Set default exchange bounds from media composition
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
    
    # Set up logistic parameters and time scale
    params = {mi: params[mn] for mn, mi in metmap.items() if mn in params}
    ts_array = np.linspace(0, t_max, num=t_max*resolution+1)
    timecourse = list(ts_array)
    for mi, param_set in params.items():
        param_set.eval(ts_array)
    for mi in "valL", "isobuta":
        params[mi].tshift(params["isocap"].halfmax())
    for mi in "ileL", "2mbut":
        params[mi].tshift(params["proL"].halfmax())
    
    # Initialize data structure for tracked metabolites
    propdata = dict()
    propmets = []
    for mi in tracked_mets:
        if model.metabolites.has_id(mi + "_c"):
            propdata[mi] = {
                'rxns_in': set([]),
                'rxns_out': set([]),
                'data_in': [],
                'data_out': [],
            }
            propmets.append(model.metabolites.get_by_id(mi + "_c"))
            if fva_run:
                propdata[mi]['data_in_lb'] = []
                propdata[mi]['data_in_ub'] = []
                propdata[mi]['data_out_lb'] = []
                propdata[mi]['data_out_ub'] = []
        else:
            print(f"'{mi}_c' is not a valid metabolite.")

    # Initialize results dataframes
    mnames = [model.metabolites.get_by_id(mn + '_c').name for mn in params]
    results = [pd.DataFrame(0., index=timecourse, columns=mnames) for _ in range(6)]
    rnames = [model.reactions.get_by_id(ri).name for ri in tracked_rxns]
    rxnflux = pd.DataFrame(0., index=timecourse, columns=rnames)
    if fva_run:
        rxnf_ub = pd.DataFrame(0., index=timecourse, columns=rnames)
        rxnf_lb = pd.DataFrame(0., index=timecourse, columns=rnames)
        allrxns = [rxn.id for rxn in model.reactions]
        fullflux_lb = pd.DataFrame(0., index=timecourse, columns=allrxns)
        fullflux_ub = pd.DataFrame(0., index=timecourse, columns=allrxns)

    print('dFBA log: simulation output begin.')
    # Simulate static solutions over timecourse
    for i, t in enumerate(timecourse):
        # Update exchange flux bounds for each constrained metabolite
        for j, mi in enumerate(params):
            ## Calculate and set exchange constraints ##
            fbounds = update_uptake_bounds(model, t, mi, params[mi])
            mn = model.metabolites.get_by_id(mi + '_c').name

            ## Record exchange constraints ##
            if not (obj == "ATP_sink" and mi == 'glc'):
                for k, df in enumerate(results):
                    df.at[t, mn] = fbounds[k]

            ## Calculate glucose uptake constraint ##
            if obj == "ATP_sink" and mi in ['ac', 'eto', 'but', 'alaL']:
                for k, df in enumerate(results):
                    if mi == 'but':
                        df.at[t, 'beta-D-glucose'] += fbounds[k]
                    else:
                        df.at[t, 'beta-D-glucose'] += 0.5*fbounds[k]

        ## Set glucose uptake constraint ##
        if obj == "ATP_sink":
            set_bounds(model, "Ex_glc", lower=results[1].at[t, 'beta-D-glucose'],
                       upper=results[2].at[t, 'beta-D-glucose'], update=True)

        ## Add natural abundance acetate allowance to secretion constraint ##
        lower, upper = model.reactions.Sec_ac.bounds
        set_bounds(model, "Sec_ac", lower=lower+results[1].at[t, 'glycine'],
                   upper=upper+results[2].at[t, 'glycine'], update=True)

        ## Get flux solution(s) and populate arrays ##
        sol = fba_method(model)
        if fva_run:
            sol_v = cb.flux_analysis.flux_variability_analysis(
                model,
                fraction_of_optimum=0.995,
                # loopless=True,
            )
            fullflux_lb.loc[t, :] = sol_v['minimum']
            fullflux_ub.loc[t, :] = sol_v['maximum']
        if i % 10 == 0:
            print(f'dFBA log: Time = {t}  (cycle {i+1}) \tFBA solution: ' \
                  f'{sol.fluxes[obj]}')
        if sol.fluxes[obj] < 0.0001:
            print(f'dFBA log: infeasible solution on cycle {i}.')
        for ri in tracked_rxns:
            rn = model.reactions.get_by_id(ri).name
            rxnflux.at[t, rn] = sol.fluxes[ri]
            if fva_run:
                rxnf_ub.at[t, rn] = sol_v.at[ri, 'maximum']
                rxnf_lb.at[t, rn] = sol_v.at[ri, 'minimum']

        ## Display incremental solutions ##
        if t in [0, 6, 8, 12, 24, 36, 48]:
            print("time = " + str(i))
            print(model.summary(solution=sol))
            for met in propmets:
                print(met.summary(solution=sol))

        ## Record flux contributions for tracked metabolites ##
        for met in propmets:
            mi = met.id[:-2] # metabolite ID, without compartment label
            flux_dic = {}   # container for metabolite flux data
            # Collect contributions for metabolite influx and outflux
            for dr in ('in', 'out'):
                # fluxes in direction <dr>
                flux_dic[f"data_{dr}"] = dict()
                if fva_run:
                    flux_dic[f"data_{dr}_lb"] = dict()
                    flux_dic[f"data_{dr}_ub"] = dict()
            # Populate, considering all scenarios where met is produced/consumed
            for rxn in met.reactions:
                # Record outflux if meets threshold
                if (met in rxn.reactants and sol.fluxes[rxn.id] > 1E-5) \
                   or (met in rxn.products and sol.fluxes[rxn.id] < -1E-5):
                    flux_dic[f"data_out"][rxn.id] = rxn.metabolites[met] \
                                                   * sol.fluxes[rxn.id]
                    if fva_run:
                        flux_dic[f"data_out_lb"][rxn.id] = rxn.metabolites[met] \
                                                          * sol_v.at[rxn.id, 'minimum']
                        flux_dic[f"data_out_ub"][rxn.id] = rxn.metabolites[met] \
                                                          * sol_v.at[rxn.id, 'maximum']
                # Record influx if meets threshold
                elif (met in rxn.products and sol.fluxes[rxn.id] > 1E-5) \
                     or (met in rxn.reactants and sol.fluxes[rxn.id] < -1E-5):
                    flux_dic[f"data_in"][rxn.id] = rxn.metabolites[met] * \
                                                    sol.fluxes[rxn.id]
                    if fva_run:
                        flux_dic[f"data_in_lb"][rxn.id] = rxn.metabolites[met] \
                                                           * sol_v.at[rxn.id, 'minimum']
                        flux_dic[f"data_in_ub"][rxn.id] = rxn.metabolites[met] \
                                                           * sol_v.at[rxn.id, 'maximum']
            for dr in ('in', 'out'):
                propdata[mi][f'rxns_{dr}'].update(flux_dic[f"data_{dr}"].keys())
                propdata[mi][f'data_{dr}'].append(flux_dic[f'data_{dr}'])
                if fva_run:
                    propdata[mi][f'rxns_{dr}'].update(flux_dic[f"data_{dr}_lb"].keys())
                    propdata[mi][f'rxns_{dr}'].update(flux_dic[f"data_{dr}_ub"].keys())
                    propdata[mi][f'data_{dr}_lb'].append(flux_dic[f'data_{dr}_lb'])
                    propdata[mi][f'data_{dr}_ub'].append(flux_dic[f'data_{dr}_ub'])

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
    print(f'Complete after {i} cycles ({t} hours). Final flux: '\
          f' {sol.fluxes["Ex_biomass"]}')
    areaplot(results[0], results[1], results[2], ylabel='flux (mM/h)')
    plt.show()
    areaplot(results[3], results[4], results[5], ylabel='normalized signal')
    plt.show()
    niceplot(rxnflux)
    plt.show()

    ## Write run results ##
    if not dry_run:
        writer = pd.ExcelWriter(f'{SCDIR}/../../data/fluxes.xlsx', engine='openpyxl')
        snames = ('data', 'lbdev', 'ubdev', 'signal', 'slbdev', 'subdev')
        # Write estimated concentrations and fluxes
        for i, df in enumerate(results):
            df.to_excel(writer, sheet_name=snames[i], engine='openpyxl')
        # Write tracked reaction fluxes
        rxnflux.to_excel(writer, sheet_name='fluxes', engine='openpyxl')
        if fva_run:
            rxnf_ub.to_excel(writer, sheet_name='fluxub', engine='openpyxl')
            rxnf_lb.to_excel(writer, sheet_name='fluxlb', engine='openpyxl')
            # Write FVA bounds for ALL reactions to tab-delimited text
            fullflux_lb.to_csv(
                f'{SCDIR}/../../data/fullfluxlb.txt',
                sep='\t',
                index_label='Time (h)'
            )
            fullflux_ub.to_csv(
                f'{SCDIR}/../../data/fullfluxub.txt',
                sep='\t',
                index_label='Time (h)'
            )
        writer.close()
        # Open workbook to record flux fractions for tracked metabolites
        writer = pd.ExcelWriter(f'{SCDIR}/../../data/met_fluxes.xlsx', engine='openpyxl')
        wb=writer.book
        ws=wb.create_sheet("MetFluxes")
        writer.sheets["MetFluxes"] = ws
        startcol = 0

    # Record reaction flux fractions for tracked metabolites
    for met in propmets:
        mi = met.id[:-2]
        print(f"Flux slices for metabolite {met.name} ({met.id}).")
        flux_dfs = {}
        rnames = {}
        rids = {}
        names = ["Fluxes", "Flux Upper Bound", "Fluxes Lower Bound"]
        for dr in ('in', 'out'):
            rids[dr] = list(propdata[mi][f'rxns_{dr}'])
            rnames[dr] = [model.reactions.get_by_id(rid).name for rid in rids[dr]]
            labels = [f'data_{dr}']
            if fva_run:
                labels.extend([f'data_{dr}_lb', f'data_{dr}_ub'])
            for i, lab in enumerate(labels):
                # Convert flux fraction records to dataframes
                flux_dfs[lab] = pd.DataFrame.from_records(
                    propdata[mi][lab], 
                    index=timecourse
                ).fillna(0.)
                # Print result
                print(names[i])
                print(f"\t{dr}flux percentages:")
                v = flux_dfs[lab]
                totalflux = v.to_numpy().sum()
                cuts = v.sum(axis=0)
                for i, rxn in enumerate(rids[dr]):
                    print("\t", rxn, rnames[dr][i], f"{cuts.at[rxn]/totalflux*100:.2f}")
                print(f"\tAbsolute {dr}flux:")
                for i, rxn in enumerate(rids[dr]):
                    print("\t", rxn, rnames[dr][i], f"{cuts.at[rxn]}")
            # Write flux fraction results to Excel (met_fluxes.xlsx)
            if not dry_run:
                df = flux_dfs[f'data_{dr}']
                endcol = startcol + df.shape[1]
                if startcol == 0:
                    write_index = True
                    endcol += 1
                    ws.cell(row=1, column=2, value="Time")
                else:
                    write_index = False
                df.to_excel(writer, sheet_name="MetFluxes", startrow=2,
                            startcol=startcol, index=write_index)
                c = ws.cell(row=1, column=startcol+1, value=f"{met.name} {dr}flux")
                c.font = BOLD
                try:
                    ws.merge_cells(start_row=1, start_column=startcol+1, end_row=1,
                                   end_column=endcol)
                except ValueError:
                    pass
                colnames = list(df.columns)
                for i, j in enumerate(range(endcol-df.shape[1], endcol)):
                    rxname = model.reactions.get_by_id(colnames[i]).name
                    c = ws.cell(row=2, column=j+1, value=rxname)
                    c.font = BOLD
                startcol = endcol
    if not dry_run:
        writer.close()

def handle_args(args):
    """Parse command line call and run dFBA function."""
    if len(args) < 1:
        print("Please provide a FBA method to use. Options: " +
              ", ".join(METHODS))
        return
    elif len(args) < 2:
        print("Please provide a FBA method and paths to at least one dataset.")
        return
    method = args.pop(0)
    fva_run = False
    if method == 'fva':
        sim_fun = METHODS['fba']
        fva_run = True
    elif method in METHODS:
        sim_fun = METHODS[method]
    else:
        print(f"Unsupported FBA method {method}. Please select one of "
                + ", ".join(METHODS))
        return
    dirs = [fp.rstrip('/') for fp in args]
    params = defaultdict(LogisticSet)
    for i, tscale in enumerate(synchronizers(dirs, stretch=stretch, \
            plot=False)):
        filepath = dirs[i]
        curves, errors = fit_trajectories(filepath, tscale=tscale, plot=True)
        for met, pset in curves.items(): # update LogisticSet of metabolite
            params[met].add_curve(pset, errors[met])
    t_max = dFBA_kwargs['t_max']
    t_num = dFBA_kwargs['t_max']*dFBA_kwargs['resolution']
    print("=========")
    print("dFBA run average curveshapes.")
    print("---------")
    for met, curveset in params.items():
        print(met)
        param_names = ["L", "k", "x0"]
        if met in ["Glucose", "Leucine", "Proline"]:
            param_names.append("C")
        curveset.display_avg_coeffs(param_names)
    areaplot2(np.linspace(0, t_max, num=t_num), params)
    dfba_main(
        params, tracked_reactions, fba_method=sim_fun, fva_run=fva_run,
        **dFBA_kwargs, tracked_mets=tracked_metabolites
    )

if __name__ == "__main__":
    """Handle command line call"""
    handle_args(sys.argv[1:])
    
