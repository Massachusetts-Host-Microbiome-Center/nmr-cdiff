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
import argparse
import datetime
import json
import os

import cobra as cb
from matplotlib import pyplot as plt
import matplotlib.colors
import numpy as np
import openpyxl as xl
import pandas as pd

from curveshapes import Metabolite
from standard_curves import compute_standard_curves
from synchronize import synchronizers
from trajectories import fit_trajectories
from get_color import get_cmap
import configparser



SCDIR = os.path.dirname(os.path.abspath(__file__))   # location of script
BOLD = xl.styles.Font(bold=True)
#MODEL_PATH = f'{SCDIR}/../data/icdf843.json'
METHODS = {
    'fba': lambda m: m.optimize(),
    'fva': cb.flux_analysis.flux_variability_analysis,
    'loopless': cb.flux_analysis.loopless.loopless_solution,
    'pfba': cb.flux_analysis.pfba,
}

class MetaboliteCollection:
    def __init__(self):
        self.name_map = {}
        self.id_map = {}

    def new(self, met_id, met_name, substrate_name, scale):
        if (met_name not in self.name_map) and (met_id not in self.id_map):
            met = Metabolite(met_id, met_name)
            self.name_map[met.name] = met
            self.id_map[met.id] = met
        else:
            met = self.get_by_id(met_id)
        met.set_substrate_scale(substrate_name, scale)

    def get(self, met_id):
        return self.get_by_id(met_id)

    def get_by_id(self, met_id):
        if met_id in self.id_map:
            return self.id_map[met_id]
        else:
            raise ValueError(f"Collection does not contain metabolite {met_id}.")
        
    def has_id(self, met_id):
        return met_id in self.id_map
    
    def remove_met(self, met_id):
        if self.has_id(met_id):
            self.pop(met_id)
 
    def get_by_name(self, name):
        if name in self.name_map:
            return self.name_map[name]
        else:
            raise ValueError(f"Collection does not contain metabolite {name}.")
        
    def get_values(self):
        return self.id_map.values()
    
    def get_ids(self):
        return self.id_map.keys()

    def get_items(self):
        return self.id_map.items()
    
    def pop(self, mid):
        item = self.id_map.pop(mid)
        self.name_map.pop(item.name)
        return item

class Substrate():
    def __init__(self, jobj):
        required_fields = ["name", "label", "model_id", "concentration", "curve", "standard_peaks", "standard_concentrations",
                           "experiments", "products", "reactions_to_constrain", "normalize_percent"]

        missing_keys = [key for key in required_fields if key not in jobj]
        if missing_keys:
            error_message = f"The following keys are missing: {missing_keys}"
            raise ValueError(error_message)


        self.name = jobj['name']
        self.label = jobj['label']
        self.model_id = jobj['model_id']
        self.concentration = jobj['concentration']
        self.curve = jobj['curve']
        self.experiments = jobj['experiments']
        self.products = jobj['products']
        self.reactions_to_constrain = jobj['reactions_to_constrain']
        self.fname_standard_peaks = jobj['standard_peaks']
        self.fname_standard_concentrations = jobj['standard_concentrations']
        self.plot_ticks = jobj['plot_ticks']
        self.cx = jobj['concentration']
        self.from_standards = False
        self.normalize_percent = jobj['normalize_percent']
        print(self.fname_standard_peaks)
        if self.fname_standard_peaks != "None":
            self.from_standards = True

        
        self.product_ids = {k: v["model_id"] for k, v in self.products.items()}
        self.product_curves = {k: v["curve"] for k, v in self.products.items()}

        if self.from_standards:
            print("Reading standards")
            fpath = parse_filepath(self.fname_standard_peaks)
            concentrations = parse_filepath(self.fname_standard_concentrations)
            
            num_carbons = {}
            for product_key, product in self.products.items():
                num_carbons[product_key] = product['num_carbons']
            sheetname = "area"
            if jobj['normalize_percent']:
                sheetname = "area_percent"
            self.product_scale = compute_standard_curves(self.fname_standard_peaks, self.fname_standard_concentrations,
            num_carbons, substrate=self.name, sheetname=sheetname)
            print("Product scale")
            print(self.product_scale)
        else:
            try:
                self.product_scale = {k: v['scale'] for k, v in self.products.items()}
            except KeyError:
                print("Expected field \"scale\" for all products of metabolite " \
                      + f"{self.name} where path to standards file was not proveded, " \
                      + "but at least one \"scale\" field was missing.")
                print("Please provide either scale fields or a standards field.")
                raise



    def metmap(self):
        met_dict = {self.name: self.model_id}
        met_dict.update({k: v for k, v in self.product_ids.items()})
        return met_dict




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
    for met, ser in df.items():
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
    plt.legend()
    plt.show()



def areaplot2(t, substrates, met_collect: MetaboliteCollection):
    """Plot lines with shaded confidence interval.

    Parameters:
    df -- dataframe of vectors corresponding to optimal values
    dl -- dataframe of vectors corresponding to lower bounds
    du -- dataframe of vectors corresponding to upper bounds
    t_max -- maximum value of x-axis
    ylabel -- y-axis label
    """

    tmax = np.max(t)
    cmap = get_cmap()
    figmap = {}
    for substrate in substrates:
        figmap[substrate.name] = plt.figure(figsize=(2.5, 2), constrained_layout=True)
        ax = figmap[substrate.name].add_subplot(111)
        afont = {'fontname': 'Arial', 'size': 7}
        ax.set_xlabel("Time (h)", **afont)
        ax.set_ylabel("Estimated Concentration (mM)", **afont)
        ax.set_xlim((0, tmax))
        if substrate.plot_ticks is not None:
            start = substrate.plot_ticks["start"]
            stop = substrate.plot_ticks["stop"]
            major = substrate.plot_ticks["major_step"]
            minor = substrate.plot_ticks["minor_step"]
            ax.set_yticks(np.arange(start, stop+major, major))
            ax.set_yticks(np.arange(start, stop+minor, minor), minor=True)
        ax.set_xticks([0, 12, 24, tmax])
        ax.set_xticks(list(range(tmax)), minor=True)
        ax.set_xticklabels(ax.get_xticks(), **afont)
        ax.set_yticklabels(ax.get_yticks(), **afont)
        #ax.xaxis.set_tick_met_collect(width=0.5)
        #ax.yaxis.set_tick_met_collect(width=0.5)
        plt.setp(ax.spines.values(), linewidth=0.5)
    
    print("PLOTTING HERE: ")
    for cpdset in met_collect.get_values():
        for substrate_name, curveset in cpdset.logistic_sets.items():
            print(cpdset.name)
            print(curveset)
            print(curveset.curves)
            f = figmap[substrate_name]
            ax = f.get_axes()[0]
            try:
                ser = curveset.get_sol(t)[0]
            except Exception as e:
                print("ERROR: Check your environment is loaded correctly. See installing_nmr_processing_environment.md")
                quit()
   
            _, _, lb, ub = curveset.get_bounds(t)
            color = cmap[cpdset.name]
            ax.plot(t, ser, '-', label=cpdset.name, lw=2, c=color)
            color = color + (0.2,)
            ax.fill_between(t, lb, ub, color=color)
            ax.legend(bbox_to_anchor = (1.01, 1.01))
    
    plt.show()
    return(f)
    # for sub, fig in figmap.items():
    #     plt.figure(fig.number)
    #     plt.show()

def load_model(modelfile, objective_list):
    """Load model and set constraints."""
    model = cb.io.load_json_model(modelfile)
    objective_dict = {}
    for objective in objective_list:
        objective_reaction = model.reactions.get_by_id(objective)
        objective_dict[objective_reaction] = 1 / len(objective_list) # assumes all equally important
    model.objective = objective_dict
    model.reactions.get_by_id(objective).upper_bound=1000

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
    model.reactions.Ex_cysL.upper_bound = 1000
    model.solver = 'glpk'
    return model


def manually_adjust_leucine_metabolism(met_collect):

    # shift the halfmax of isocaproate curve to match the halfmax of the Leucine curve (leucine is directly measured)
    # shift the halfmax of the valine curve to the same.
    # if there is no isocaproate (had KO), then just shift the valine curve
    if met_collect.has_id("leuL") and met_collect.has_id("valL") and met_collect.has_id("isobuta"):
        for mi in "valL", "isobuta": # Remove these 4 lines if we get Val and Ile runs
            try:
                met_collect.get_by_id(mi).tshift("Leucine", met_collect.get_by_id("isocap").avg_x0("Leucine"))
                met_collect.get_by_id(mi).tshift("Leucine", met_collect.get_by_id("ival").avg_x0("Leucine"))
            except:
                met_collect.get_by_id(mi).tshift("Leucine", met_collect.get_by_id("ival").avg_x0("Leucine"))
    else:
        met_collect.remove_met("valL")
        met_collect.remove_met("isobuta")
    if met_collect.has_id("proL") and met_collect.has_id("ileL") and met_collect.has_id("2mbut"):
        for mi in "ileL", "2mbut":
            met_collect.get_by_id(mi).tshift("Leucine", met_collect.get_by_id("proL").avg_x0("Proline"))
    else:
        met_collect.remove_met("ileL")
        met_collect.remove_met("2mbut")
    return(met_collect)



def initialize_result_storage(model, met_collect, timecourse, tracked_reactions, tracked_metabolites):
    # Initialize data structure for tracked metabolites
    print("Initialize data structure for tracked metabolites")
    print("Planning to track: ", tracked_metabolites)
    propdata = dict()
    propmets = []
    for mi in tracked_metabolites:
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
            print(f"'Tracked metabolite {mi}_c' is not a valid metabolite.")

    # Initialize results dataframes for trakced reactions
    mnames = [model.metabolites.get_by_id(mi + '_c').name for mi in met_collect.get_ids()]
    results = [pd.DataFrame(0., index=timecourse, columns=mnames) for _ in range(6)]
    tracked_reactions_names = [model.reactions.get_by_id(ri).name for ri in tracked_reactions]
    rxnflux_tracked = pd.DataFrame(0., index=timecourse, columns=tracked_reactions_names)
    
    # Initialize results dataframes for all reactions
    rnames_all = [model.reactions.get_by_id(ri.id).name for ri in model.reactions]
    metnames_all = [met.name for met in model.metabolites]
    rxnflux_all = pd.DataFrame(0., index=timecourse, columns=rnames_all)
    allrxns = [rxn.id for rxn in model.reactions]

    # Initialize results dataframes for upper/lower bound of reactions (used for fva)
    rxnf_ub = pd.DataFrame(0., index=timecourse, columns=tracked_reactions_names)
    rxnf_lb = pd.DataFrame(0., index=timecourse, columns=tracked_reactions_names)
    fullflux_lb = pd.DataFrame(0., index=timecourse, columns=allrxns)
    fullflux_ub = pd.DataFrame(0., index=timecourse, columns=allrxns)

    return(propdata, propmets, results, rxnflux_tracked, rxnflux_all, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb)



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



def update_uptake_bounds(model, t, met, curveset, update=True):
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
    # Calculate logistic estimates of metabolite concentration at time t
    signal, serr, exch, exerr = curveset.get_sol(t)
    exch_l, exch_u, signal_l, signal_u = curveset.get_bounds(t)

    # Reverse direction and set bounds for secretion reactions
    if met in ['glc', 'proL', 'valL', 'ileL', 'thrL']:
        rid = 'Ex_' + met
        exch, exch_l, exch_u = reverse_flux(exch, exch_l, exch_u)
        set_bounds(model, rid, lower=exch_l, upper=exch_u, update=update)

    if met in ['leuL']:
        set_bounds(model, rid, update=update) # leave leucine unbounded

    # Ignore data from 1H spectra
    elif met == 'acoa':
        pass
    # Update bounds for products
    elif met in ['2abut', 'ppa']:
        rid = 'Sec_' + met
        set_bounds(model, rid, lower=exch_l, upper=exch_u, update=update)
    else:
        rid = 'Sec_' + met
        lb = max(exch_l, 0) # do not allow reverse flux
        ub = max(exch_u, 0)
        set_bounds(model, rid, lower=lb, upper=ub, update=update)
    # Update Wood-Ljungdahl Pathway bounds with butyrate
    if met == 'but':
        wlp_l, wlp_u, _, _ = curveset.get_bounds(t, substrate="Glucose")
        set_bounds(model, 'ID_326', lower=0, upper=max(wlp_u, 0), update=update)
    # Allow natural abundance acetate from cysteine
    #if met == 'ac':
    #    cys_l, cys_u, _, _ = curveset.get_bounds(t, substrate="Acetate13C")
    #    set_bounds(model, 'Ex_cysL', lower=max(cys_l, 0), upper=max(cys_u, 0), update=update)


    return exch, exch_l, exch_u, signal, signal_l, signal_u



def manually_constrain_reaction_to_substrate(t, substrate, rid, curveset, model):
    # in specific cases, we know that the activity of a particular reaction need be constrained by the transport/availability of one of the substrates
    # for example, part of the butyrate biosynthesis pathway ID_325 enoyl-CoA hydratase is constrained by glucose uptake
    # and butyrate biosynthesis pathway reaction 2HBD is instead constrained by Threonine
    if substrate in curveset.logistic_sets:
        lower_bound, upper_bound, _, _ = curveset.get_bounds(t, substrate=substrate)
        set_bounds(model, rid, lower=lower_bound, upper=upper_bound, update=True)
    return(model)


def get_flux_contributions_to_metabolites(propdata, propmets, sol):
    ## Record flux contributions for tracked metabolites ##
    # propmets = list of cobra metabolite objects
    print("Prop mets: ", propmets)
    for met in propmets:
        mi = met.id.replace("_c", "").replace("_e", "")
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
    return(propdata)





def dfba_main(met_collect: MetaboliteCollection, model_file, objective_function, fba_method, substrates, seed,
              fva_run=False, tracked_reactions=[], tracked_metabolites=[], tmin_hours=0,
              tmax_hours=48, solutions_per_hour=1, dry_run=False, output_folder = "../data/", reactions_to_delete = []):
    
    print("Tracked reactions: ", tracked_reactions)
    print(met_collect.get_items())

    """Main function to compute dFBA solutions.
    Computes successive static FBA solutions and plots the estimated metabolite
    concentrations, uptake rates, and tracked reaction fluxes.

    Parameters:
    met_collect -- dictionary mapping metabolite names to LogisticSet objects
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

    np.random.seed(seed)
    print(f"""dFBA log: Begin dFBA analysis with sheet, endpoint
          {tmax_hours} hours, and resolution {solutions_per_hour}.""")
    print(f"Using method {fba_method}.")
    if dry_run:
        print("Dry run, will not write results.")
    # Load metabolic model and logistic fit specs #
    print('dFBA log: loading model and specsheet...')
    model = load_model(parse_filepath(model_file), objective_function)


    # Remove reactions to delete
    model.remove_reactions(reactions_to_delete)
    
    # 1. Evaluate the learned logistic function per metabolite for every timepoint
    nsol = int(round((tmax_hours - tmin_hours)*solutions_per_hour + 1, 0))
    ts_array = np.linspace(tmin_hours, tmax_hours, num=nsol)
    timecourse = list(ts_array)
    for mi, param_set in met_collect.get_items():
        param_set.eval(ts_array)

    # Special-case adjustments for certain metabolites
    #met_collect.remove_met("5apn") #Aidan, why did yuo remove this proline metabolite?

    # 2. Manually adjust leucine metabolism
    met_collect = manually_adjust_leucine_metabolism(met_collect)
    

    # 3. initialize results storage
    propdata, propmets, results, rxnflux_tracked, rxnflux_all, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb = initialize_result_storage(model, met_collect, timecourse, tracked_reactions, tracked_metabolites)
    allrxns = [rxn.id for rxn in model.reactions]

    print('dFBA log: simulation output begin.')
    # Simulate static solutions over timecourse

    for i, t in enumerate(timecourse):

        # 1. Update exchange flux bounds for each constrained metabolite
        for j, (mi, curveset) in enumerate(met_collect.get_items()):
            
            # Use NMR data to constrain boundary reactions. This will inherently limit the flux that may be assigned to transport reactions involving this metabolite
            fbounds = update_uptake_bounds(model, t, mi, curveset)

            ## Record exchange constraints
            mn = model.metabolites.get_by_id(mi + '_c').name
            if not (objective_function == "ATP_sink" and mi == 'glc'):
                for k, df in enumerate(results):
                    df.at[t, mn] = fbounds[k]

            # manually constrain butyrate to be produced by only glucose, or only threonine if threonine is measured by NMR
            if mi == 'but':
                model = manually_constrain_reaction_to_substrate(t, "Glucose", "ID_325", curveset, model)
                model = manually_constrain_reaction_to_substrate(t, "Threonine", "2HBD", curveset, model)


        #2/ Get flux solution(s) and populate arrays ##
        print("Starting solving")
        sol = fba_method(model)

        print("Objective function: ", sol.fluxes[objective_function])
        if fva_run:
            print("FVA")
            sol_v = cb.flux_analysis.flux_variability_analysis(
                model,
                fraction_of_optimum=0.995,
                loopless=False
                # loopless=True
            )
            fullflux_lb.loc[t, :] = sol_v['minimum']
            fullflux_ub.loc[t, :] = sol_v['maximum']

        # 2.5 print out results for verbose logging
        if i % 10 == 0:
            print(f'dFBA log: Time = {t}  (cycle {i+1}) \tFBA solution: ' \
                  f'{sol.fluxes[objective_function]}')

        if np.sum(sol.fluxes[objective_function]) < 0.0001:
            print(f'dFBA log: infeasible solution on cycle {i}.')

        # 3/ Record reaciton fluxes at time t
        for ri in tracked_reactions:
            rn = model.reactions.get_by_id(ri).name
            rxnflux_tracked.at[t, rn] = sol.fluxes[ri]
            if fva_run:
                rxnf_ub.at[t, rn] = sol_v.at[ri, 'maximum']
                rxnf_lb.at[t, rn] = sol_v.at[ri, 'minimum']

        for ri in allrxns:
            rn = model.reactions.get_by_id(ri).name
            rxnflux_all.at[t, rn] = sol.fluxes[ri]


        ## 4/ Print out log for tracking purposes ##
        if t in [0, 6, 8, 12, 21, 24, 36, 48]:
            print("Display incremental solutions")
            print("time = " + str(i))
            print(model.summary(solution=sol))
            for met in propmets:
                print(met.summary(solution=sol))


        # FVA takes longer, print at end of each simulation
        if fva_run:
            now = datetime.datetime.now()
            print(f"t = {t:.2f} simulation complete ({now:%c})")

        ## 5/ Track relative contributions to metabolites in "tracked_metabolites" list from configuration file
        propdata = get_flux_contributions_to_metabolites(propdata, propmets, sol)    

    return(propdata, propmets, results, rxnflux_tracked, rxnflux_all, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb, model)





def read_and_validate_cfg(cfg_filepath):
    parsed_path = parse_filepath(cfg_filepath)
    with open(parsed_path, "r") as rf:
        cfg = json.loads(rf.read())
    required_fields = ["method", "model_file", "objective_function", "nmr_substrates"]
    for field in required_fields:
        try:
            val = cfg[field]
        except KeyError:
            print(f"JSON config missing required field {field}.")
    return cfg

def parse_filepath(fpath):
    if not fpath.startswith('/'):
        return f"{SCDIR}/{fpath}"
    return fpath


def write_out_fluxes(output_folder, results, rxnflux_all, rxnflux_tracked, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb):
    writer = pd.ExcelWriter(f"{output_folder}/fluxes.xlsx", engine='openpyxl')
    snames = ('data', 'lbdev', 'ubdev', 'signal', 'slbdev', 'subdev')
    # Write estimated concentrations and fluxes
    for i, df in enumerate(results):
        df.to_excel(writer, sheet_name=snames[i], engine='openpyxl')
    # Write tracked reaction fluxes
    rxnflux_tracked.to_excel(writer, sheet_name='fluxes', engine='openpyxl')
    if fva_run:
        rxnf_ub.to_excel(writer, sheet_name='fluxub', engine='openpyxl')
        rxnf_lb.to_excel(writer, sheet_name='fluxlb', engine='openpyxl')
        # Write FVA bounds for ALL reactions to tab-delimited text
        fullflux_lb.to_csv(
            f'{output_folder}/fullfluxlb.txt',
            sep='\t',
            index_label='Time (h)'
        )
        fullflux_ub.to_csv(
            f'{output_folder}/fullfluxub.txt',
            sep='\t',
            index_label='Time (h)'
        )
    rxnflux_all.to_excel(writer, sheet_name = "allfluxes", engine = 'openpyxl')
    writer.close()


def write_flux_contributions_to_metabolites(output_folder, propdata, propmets, model, timecourse):
    # Open workbook to record flux fractions for tracked metabolites
    writer = pd.ExcelWriter(f'{output_folder}/met_fluxes.xlsx', engine='openpyxl')
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

    writer.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--params_file', metavar='CONFIG', help='path to dFBA config file')
    #parser.add_argument('--output_folder', default = None, help = "folder in which to output all results, to maintain different tests")
    parser.add_argument('--trajectory_only', action = 'store_true', help = "True if you only want to plot the measured metabolite trajectories and not run FBA")

    args = parser.parse_args()

    # 1. Retrieve info from config file
    ### retrieve parameters
    cfg_path = args.params_file
    trajectory_only = args.trajectory_only

    ### declare a metabolite collection object
    met_collect = MetaboliteCollection()

    ### get info from config file
    cfg = read_and_validate_cfg(cfg_path)
    reactions_to_delete = [i for i in cfg.pop('reactions_to_delete')]
    print("Removing reactions: ", reactions_to_delete)
    tracked_reactions = [i for i in cfg.pop('tracked_reactions')]
    tracked_metabolites = [i for i in cfg.pop('tracked_metabolites')]



    substrates = [Substrate(s) for s in cfg['nmr_substrates']]

    isotope = cfg['isotope']
    
    
    tmax = cfg['tmax_hours']
    output_folder = cfg['output_folder']
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # initialize the metabolite logistic function holding object
    for substrate in substrates:
        met_collect.new(substrate.model_id, substrate.name, substrate.name, 1.)
        for metname, metdata in substrate.products.items():
            met_collect.new(metdata["model_id"], metname, substrate.name, substrate.product_scale[metname])


    method = cfg["method"]
    plot = cfg['plot']
    seed = cfg['seed']


    # 2. Set up objects
    ### set up experiment objects
    all_experiments = {}
    for substrate in substrates:
        all_experiments.update({parse_filepath(exp): substrate for exp in substrate.experiments})

    ### set optimization method
    if method == 'fva':
        fba_method = METHODS['fba']
        fva_run = True
    elif method in METHODS:
        fba_method = METHODS[method]
        fva_run = False
    else:
        print(f"Unsupported FBA method {method}. Please select one of "
                + ", ".join(METHODS))


    # 3. Synchronize all runs in time using start of isocaproate formation
    sync_functions = synchronizers(
        [path for path in all_experiments], 
        stretch=cfg.pop("stretch", False), 
        plot=False
    )

    for tscale, (fname_exp, substrate) in zip(sync_functions, all_experiments.items()):
        print("FNAME_EXP: ", fname_exp)
        sheetname = "area"
        if substrate.normalize_percent:
            sheetname = "area_percent"
        time, signals, curves, errors = fit_trajectories(fname_exp, isotope, tmax, substrate, sheetname = sheetname, tscale=tscale, plot=True)


        # check that if a product is listed in the expected products portion of the json file, it has associated NMR data
        missing_keys = [key for key in substrate.products if key not in list(curves.keys())]
        if missing_keys:
            error_message = f"{missing_keys} is in the expected product list in your dfba_cfg.json file, but there is no associated NMR data in {fname_exp}_13C.xlsx"
            raise ValueError(error_message)

        df = pd.DataFrame(signals, index = time)
        df.to_csv(output_folder + "/time_norm_areas.csv", index_label = "Time")

        print("Curves: ", curves)
        for met, pset in curves.items(): # update LogisticSet of metabolite
            print("Adding a logistic curve object for: ", met)
            met_collect.get_by_name(met).add_curve(substrate.name, pset, errors[met])

    if plot:
        plt.tight_layout()
        plt.savefig(output_folder + "/trajectories.png")
        print("Outputting trajectory plot to: ", output_folder + "/trajectories.png")
        plt.show()

    # set the number of experiments (nmr runs) representing every metabolite
    for substrate in substrates:
        for met in substrate.metmap():
            met_collect.get_by_name(met).set_runcount(substrate.name, len(substrate.experiments))


    # 4. fit logistic parameters to each metabolite curve, save parameters to output_folder/
    logistic_df = []
    for met_id, curveset in met_collect.get_items():
        print(met_id)
        print("\n")
        avg_ps, avg_err = curveset.display_avg_coeffs()
        param_names = ["L", "k", "x0", "C"]
        print(len(avg_ps))
        print(len(param_names))
        logistic_df.append(pd.DataFrame({'parameter': param_names[0:len(avg_ps)], 'values': avg_ps, 'met_id': met_id, 'met_name': curveset.name}))       
    pd.concat(logistic_df).to_csv(output_folder + "/logistic_met_collect.csv")

    t_min = cfg["tmin_hours"]
    t_max = cfg["tmax_hours"]
    t_num = int(round((t_max-t_min)*cfg["solutions_per_hour"] + 1, 0))

    if plot:
        if len(substrates) > 1:
            f = areaplot2(np.linspace(t_min, t_max, num=t_num), substrates, met_collect)
            plt.tight_layout()
            plt.savefig(f, output_folder + "/trajectories.png")
            print("Outputting trajectory plot to: ", output_folder + "/trajectories.png")




    # 5. Run dfba
    objective = cfg['objective_function'][0]
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    if not trajectory_only:
        propdata, propmets, results, rxnflux_tracked, rxnflux_all, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb, model = dfba_main(met_collect, 
            objective_function = cfg.pop('objective_function'),
            model_file = cfg.pop('model_file'),
            fva_run=fva_run, fba_method=fba_method,
         substrates=substrates, seed = seed, output_folder = output_folder,
          reactions_to_delete = reactions_to_delete,
           tracked_reactions = tracked_reactions, tracked_metabolites = tracked_metabolites)   


        # 6. Write out results
        ## Adjust flux solution directionality in output for visualization ##
        ## TO BE REMOVED AFTER CHECKING WITH UNIDIRECTIONAL MODEL ##
        reverse_ids = ["ID_383", "ID_336", "ID_391", "HydEB"]
        reverse_nms = [model.reactions.get_by_id(rid).name for rid in reverse_ids]
        for rname in reverse_nms:
            if rname in rxnflux_tracked:
                rxnflux_tracked[rname] *= -1
                if fva_run:
                    rxnf_ub[rname] *= -1
                    rxnf_lb[rname] *= -1

        ## Plot run results ##
        #print(f'Complete after {i} cycles ({t} hours). Final flux: '\
        #      f' {sol.fluxes["Ex_biomass"]}')
        #areaplot(results[0], results[1], results[2], ylabel='flux (mM/h)')
        #areaplot(results[3], results[4], results[5], ylabel='normalized signal')
        #niceplot(rxnflux_tracked)

        ## Write run results ##
        # 1. write out fluxes (Tracked and all)
        write_out_fluxes(output_folder, results, rxnflux_all, rxnflux_tracked, rxnf_ub, rxnf_lb, fullflux_ub, fullflux_lb)


        # 2. Write out flux contributions to tracked metabolites
        print(rxnflux_all.index.values)
        write_flux_contributions_to_metabolites(output_folder, propdata, propmets, model, timecourse = list(rxnflux_all.index.values))

#all_experiments: dict['nmr_data_filepath'] = Substrate()
#Substrate() class: all elements of json config file