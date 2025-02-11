import pandas as pd
import numpy as np
import modelseedpy 
import json
import re
from bioservices import KEGG
import re
from cobra.io import load_json_model, save_json_model
from cobra import Model, Reaction, Metabolite
from bioservices import KEGG
import re
from scipy.stats import ranksums
import statsmodels.stats.multitest as multi
import math
from copy import deepcopy
from cobra.flux_analysis import flux_variability_analysis
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
#import pygraphviz
import matplotlib.pyplot as plt
from cobra.medium import minimal_medium
from cobra.flux_analysis import find_essential_reactions
from collections import Counter
from scipy.stats import spearmanr
from pyvis.network import Network

from scipy.optimize import minimize
from skopt import gp_minimize
from skopt.space import Categorical
from skopt.utils import use_named_args

class Model:
    def __init__(self, model, config):
        self.model = model
        self.reaction_id_list = [i.id for i in model.reactions]
        self.metabolite_id_list = [i.id for i in model.metabolites]
        self.reaction_name_list = [i.name for i in model.reactions]
        self.metabolite_name_list = [i.name for i in model.metabolites]

        project_dir = config.get('paths', 'project_dir')
        database_path = project_dir + config.get('paths', 'modelseed_reaction_db')
        self.db = pd.read_csv(database_path)

        self.model = self.addCompartments()
        self.transport_reactions = self.getTransportReactions()
        self.medium = model.medium
        self.minimal_medium = None
        self.medium_search_space = None

        # cdiff model conversions
        compound_path1 = project_dir + "/" + config.get('paths', 'modelseed_compound_db')
        self.db_compound = get_db_compound(compound_path1)

        compound_path2 = project_dir + "/" + config.get("paths", 'modelseed_compound_db2')
        cdiff_model_file = project_dir + "/" + config.get("paths", "model_cdiff")
        self.conversions = getConversionFromCdiff(compounds_filepath = compound_path1,
                                                curated_file_path = compound_path2,
                                                cdiff_model_filepath = cdiff_model_file)
        self.boundary_to_transport = {}

        self.G = None

        self.gene_rule_df = None
        self.reachable_reaction_ids = []
        self.reachable_metabolite_ids = []



    def get_id_from_name(self, name):
        j = None
        if name in self.reaction_name_list:
            j = [i.id for i in self.model.reactions if i.name == name][0]
        if name in self.metabolite_name_list:
            j = [i.id for i in self.model.metabolites if i.name == name][0]
        return(j)


    def getFluxes(self):
        fluxes = self.model.optimize().fluxes
        fluxes = pd.DataFrame(fluxes[fluxes != 0])
        fluxes_sort = fluxes.sort_values(by = "fluxes", ascending = True)
        #tmp = self.db.loc[[i.replace("_c0", "").replace("_e0", "").replace("_z0", "") for i in fluxes_sort.index.values], :]
        #tmp = tmp.loc[[not i for i in tmp.index.duplicated(keep = "first")], :]
        #fluxes_sort['equation'] = tmp['definition'].values
        #fluxes_sort['name'] = tmp['name'].values
        #fluxes_sort['ec'] = tmp['EC_from_db'].values
        return(fluxes_sort)


    def getTransportReactions(self):
        potential_transport_reaction_ids = [i.id for i in self.model.reactions]
        keep = []
        for reaction_id in potential_transport_reaction_ids:
            reactants = self.model.reactions.get_by_id(reaction_id).reactants
            reactant_ids = [i.id for i in reactants]
            external = ["_e" in i for i in reactant_ids]
            if sum(external) > 0:
                keep.append(True)
            else:
                keep.append(False)

        transport_reactions = np.array(potential_transport_reaction_ids)[keep]

        return(transport_reactions)

    def enableTransportReaction(self, reaction_id, concentration = 30):
        self.model = self.addCompartments()
        media_dict = self.model.medium
        for reactant in self.model.reactions.get_by_id(reaction_id).reactants:
            if "_e" in reactant.id:
                new_id = "EX_" + reactant.id
                if not self.model.reactions.has_id(new_id):
                    self.model.add_boundary(self.model.metabolites.get_by_id(reactant.id), type = "exchange")
                    self.boundary_to_transport[new_id] = reaction_id
                media_dict[new_id] = concentration
        self.model.medium = media_dict


    def enableTransportReactions(self, transport_reactions = None, concentration = 30):
        transport_reactions = self.getTransportReactions()
        if transport_reactions is None:
            transport_reactions = self.transport_reactions
        for reaction_id in transport_reactions:
            self.enableTransportReaction(reaction_id, concentration)



    def addCompartments(self):
        for metabolite in self.model.metabolites:
            if "_c" in metabolite.id:
                metabolite.compartment = "c"
            if "_e" in metabolite.id:
                metabolite.compartment = "e"
        return(self.model)


    def addSink(self, product_id):
        if not self.model.reactions.has_id("SK_" + product_id):
            self.model.add_boundary(self.model.metabolites.get_by_id(product_id), type = "sink")


    def gene_rule_present(self, gene_rules):
        cdiff_genes_present = self.gene_rule_df['CD630_homolog'].values
        if isinstance(gene_rules, list):
            for rule in gene_rules:
                if rule in cdiff_genes_present:
                    return(True)
        else:
            if gene_rules in cdiff_genes_present:
                return(True)

        return(False)


    def addCdiffReaction(self, cdiff_model, cdiff_reaction_id, enable_exchange = False, check_gene_rule = True, update_reachable = True):
        #print([i.id for i in cdiff_model.reactions.get_by_id(reaction_id).metabolites.keys()])
        reaction = cdiff_model.model.reactions.get_by_id(cdiff_reaction_id)
        reaction_name = reaction.name
        gene_rules = reaction.gene_reaction_rule
        if "and" in gene_rules:
            gene_rules = gene_rules.strip().split(" and ")
        if "or" in gene_rules:
            gene_rules = gene_rules.strip().split(" or ")

        # If we're checking for gene rule, only continue forward if the gene rule is present
        continue_forward = True
        if check_gene_rule:
            if not self.gene_rule_present(gene_rules):
                #print("Gene ", gene_rules, " NOT present")
                continue_forward = False
            #else:
                #print("Gene ", gene_rules, " is present")



        reaction_return = None
        new_reaction = False
        if continue_forward:
            compound_ids = [i.id for i in reaction.metabolites.keys()]
            coefficients = [i for i in reaction.metabolites.values()]


            metabolite_set = {}
            for compound_id, coefficient in zip(compound_ids, coefficients):

                seed_id = self.conversions[compound_id]

                if seed_id not in self.db_compound.index.values:
                    self.db_compound = addDummyDatabaseEntry(seed_id, self.db_compound)

                if self.model.metabolites.has_id(seed_id):
                    met = self.model.metabolites.get_by_id(seed_id)
                else:
                    name = self.db_compound.loc[seed_id, :].values[1]
                    if name == "None" or name is None:
                        name = seed_id
                    met = Metabolite(id=seed_id, name = name)

                metabolite_set[met] = coefficient
            reaction_new = Reaction(id = cdiff_reaction_id, name = reaction_name, lower_bound=reaction.lower_bound, upper_bound=reaction.upper_bound)
            reaction_new.add_metabolites(metabolite_set)
            if not self.model.reactions.has_id(reaction_new.id):
                self.model.add_reactions([reaction_new])
                reaction_return = reaction_new
                new_reaction = True



            # If we just added a transport reaction, we need to add a corresponding exchange reaction
            # to let the model know it can get the met_e external metabolite from the media, if it is there
            cdiff_transport_reactions = cdiff_model.getTransportReactions()
            if reaction.id in cdiff_transport_reactions:
                enable_exchange = True


            if enable_exchange:
                self.model = self.addCompartments()
                self.enableTransportReaction(reaction_new.id)
                print("Enabling exchange on :", reaction_new.id)


        if reaction.id in [i.id for i in self.model.reactions]:
            reaction = self.model.reactions.get_by_id(reaction.id)
            reaction_return = reaction


        # update reachable reactions if we added a new reaction
        if update_reachable:
            transport_reactions = self.getTransportReactions()
            self.enableTransportReactions()
            reachable_reaction_ids, not_reachable_reaction_ids, reachable_metabolite_ids, not_reachable_metabolite_ids = self.findUnreachables(reachable_reaction_ids= transport_reactions)
            self.reachable_reaction_ids = reachable_reaction_ids
            print("Len of reachable reaction ids: ", len(reachable_reaction_ids), " after adding ", reaction.id)
            self.reachable_metabolite_ids = reachable_metabolite_ids




        return(reaction_return)

    def addCdiffTransportReactions(self, cdiff_model, transport_reactions_ids = None):

        if transport_reactions_ids is None:
            transport_reactions = [i  for i in cdiff_model.reactions if "Trans" in i.id or "Ex" in i.id]
            transport_reactions_ids = [i.id for i in transport_reactions]
        else:
            transport_reactions = [cdiff_model.reactions.get_by_id(i) for i in transport_reactions_ids]

        i = 0
        reactions_added = []
        for reaction_id in transport_reactions_ids:
            reaction = self.addCdiffReaction(cdiff_model, reaction_id, enable_exchange = True)
            reactions_added.append(reaction)
            i = i + 1
        return(reactions_added)

    def checkReactionObj(self, reaction_id):
        model_test = deepcopy(self)
        if "rev" in reaction_id:
            reaction_id_rev = reaction_id.replace("_rev", "")
        else:
            reaction_id_rev = reaction_id + "_rev"
        # 1. set objective
        model_test.model.objective = reaction_id
        
        # 2. add sinks for products
        products = model_test.model.reactions.get_by_id(reaction_id).products
        for product in products:
            model_test.addSink(product.id)
        
        # 3. disable the reverse to avoid infinite loop
        if model_test.model.reactions.has_id(reaction_id_rev):
            model_test.model.remove_reactions([reaction_id_rev])
        print(model_test.model.optimize().objective_value)


    def checkReactionLinked(self, reaction_id, substrate_ids, concentrations = np.arange(0, 1, 1)):
        print("Testing reaction link to substrate")
        model_test = deepcopy(self)
        if "rev" in reaction_id:
            reaction_id_rev = reaction_id.replace("_rev", "")
        else:
            reaction_id_rev = reaction_id + "_rev"
        # 1. set objective
        model_test.model.objective = reaction_id
        
        # 2. add sinks for products
        products = model_test.model.reactions.get_by_id(reaction_id).products
        for product in products:
            model_test.addSink(product.id)
        
        # 3. disable the reverse to avoid infinite loop
        if model_test.model.reactions.has_id(reaction_id_rev):
            model_test.model.remove_reactions([reaction_id_rev])
        
        # Changing leucine concentrations changes reaction id - conclusion: linked up correctly
        medium = model_test.model.medium
        obj_sequence = []
        for i in concentrations:
            for substrate_id in substrate_ids:
                medium['EX_' + substrate_id.replace("_c0", "_e0")] = i
            model_test.model.medium = medium
            val = model_test.model.optimize().objective_value
            obj_sequence.append(val)
            print("Biomass: ", val)
            #fluxes = model_test.model.optimize().fluxes
            #print(fluxes.loc[fluxes > 0])
        if spearmanr(concentrations, obj_sequence)[1] < .05:
            print("Reaction flux is correlated is substrate concentration")
        else:
            print("Reaction flux independent of substrate concentration")




    def gapfill_from_cdiff_reaction(self, cdiff_model, reaction_id, check_gene_rule):
        
        # 1. start with the unreachable reaction
        metabolites = [i.id for i in self.model.reactions.get_by_id(reaction_id).reactants]

        # 2. identify which reactants in the reaction are unreachable
        unreachable_met_ids = np.array(metabolites)[[not self.metabolite_is_reachable(i) for i in metabolites]]

        # 3. identify how template model is getting those metabolites: id reactions where product is the metabolite in question
        added_reactions = []
        met_not_reachable = [True] * len(unreachable_met_ids)
        reaction_options = {}
        for j,met in enumerate(unreachable_met_ids):
            #print(met)
            #print("unreachable metabolites: ", met, model.metabolites.get_by_id(met).name)
            reaction_options[met] = []
            for reaction in cdiff_model.model.reactions:
                product_ids = [i.id for i in reaction.products]
                product_ids_seed = [self.conversions[i] for i in product_ids if i in self.conversions]
                if met in product_ids_seed:
                    met_cdiff = [i for i in product_ids if self.conversions[i] == met]
                    print("Considering adding reaction: ", reaction.id, reaction.name, " to reach metabolite: ", met, met_cdiff)
                    reaction_options[met].append(reaction.id)

                    reaction = self.addCdiffReaction(cdiff_model, reaction.id, check_gene_rule = check_gene_rule)

                    if reaction is not None:
                        added_reactions.append(reaction.id)
                        met_not_reachable[j] = False
                        print("Metabolite : ", met, " is now reachable")

        for met, options in reaction_options.items():
            print("Reaction options for reaching: ", met)
            for option in options:
                print("\t - ", option, cdiff_model.model.reactions.get_by_id(option).name, cdiff_model.model.reactions.get_by_id(option).gene_reaction_rule)
            print("")
        return(added_reactions, unreachable_met_ids[met_not_reachable], reaction_options)


    def find_producing_reactions(self, met_id):

        met = self.model.metabolites.get_by_id(met_id)

        # Find reactions that produce these reactants
        producing_reactions = set()

        for producing_reaction in met.reactions:
            if met in producing_reaction.products:
                producing_reactions.add(producing_reaction)
        return producing_reactions


    def gapfill_from_cdiff_tree(self, cdiff_model, reaction_id_reach, cycle_lim = 3, check_gene_rule = True):



        reactions_to_add = [reaction_id_reach]

        reactions_added_to_model = []
        mets_need_to_reach = [i.id for i in self.model.reactions.get_by_id(reaction_id_reach).metabolites]
        to_break = False
        unreachable_met_ids_save = []
        while not reaction_id_reach in self.reachable_reaction_ids and len(reactions_to_add) > 0:

            print("Need to gapfill to reach reaction ids: ", reactions_to_add)

            for reaction_id in reactions_to_add:

                # add the reaction if possible, and return the metabolites that are now unreachable for the inputs of that reaction
                added_reaction_ids, unreachable_met_ids, reaction_options = self.gapfill_from_cdiff_reaction(cdiff_model, reaction_id, check_gene_rule)

                for met_id in np.array(unreachable_met_ids):
                    unreachable_met_ids_save.append(met_id)

                mets_need_to_reach = mets_need_to_reach + list(np.unique(unreachable_met_ids))
                
                reactions_to_add = []
                print("Added reaction ids: ", added_reaction_ids)
                for added_reaction_id in added_reaction_ids:
                    if added_reaction_id not in self.reachable_reaction_ids:
                        #if added_reaction_id not in reactions_added_to_model:
                        reactions_to_add.append(added_reaction_id)
                reactions_added_to_model = reactions_added_to_model + added_reaction_ids

                print(Counter(mets_need_to_reach))
                print("Compound added: ", Counter(mets_need_to_reach).most_common(1)[0][1],  ' times')
                if Counter(mets_need_to_reach).most_common(1)[0][1] > cycle_lim:
                    print("TOO MANY CYCLES")
                    to_break = True
                break;
            if to_break:
                break;
        print(Counter(mets_need_to_reach))
        print("Still unreachable mets: ", np.unique(unreachable_met_ids_save))

        return([self.model.reactions.get_by_id(i) for i in reactions_added_to_model])




    def gapfill_from_cdiff_direct(self, cdiff_model, cdiff_reaction_id, cdiff_product_id):
        cdiff_model.convertReactionsUnidirectional()

        # Add a sink for the reaction in the cdiff model
        cdiff_model.model.objective = cdiff_reaction_id
        if not cdiff_model.model.reactions.has_id("SK_" + cdiff_product_id):
            cdiff_model.model.add_boundary(cdiff_model.model.metabolites.get_by_id(cdiff_product_id), type = "sink")
        assert cdiff_model.model.optimize().objective_value > 0


        # add direct import reactions for all of the important metabolites, so that the only reactions carrying flux are those 
        # directly related to the reaction of interest
        met_ids = ['atp_c', 'adp_c', 'nad_c', 'nadh_c', 'nadph_c', 'nadp_c', 'ppi_c', 'h2o_c', 'pi_c', 'utp_c', 'ctp_c',
              'gtp_c', 'coa_c', 'fad_c', 'acoa_c', 'h_c', 'pmf_c', 'pyr_c', 'pepyr_c', 'feroxoxi_c', 'feroxred_c', 'co2_c',
              'for_c', 'h2_c', 'fmn_c', 'ribflv_c']

        medium = cdiff_model.model.medium
        for met_id in met_ids:
            if not cdiff_model.model.reactions.has_id('import_' + met_id):
                cdiff_model.addImportReaction(reaction_id = 'import_' + met_id,
                                             reaction_name = 'import_' + met_id,
                                             ex_metabolite_id = met_id.replace("_c", "_e"),
                                             in_metabolite_id = met_id,
                                             ex_metabolite_name = cdiff_model.model.metabolites.get_by_id(met_id).name,
                                             in_metabolite_name = cdiff_model.model.metabolites.get_by_id(met_id).name,
                                             enable_exchange = True,
                                             exchange_lim = 10000)

        #minimal_media_metabolites = cdiff_model.calculateMinimalMedia()
        #print(minimal_media_metabolites)
        #medium = cdiff_model.model.medium
        #for key in medium.keys():
        #    if key not in minimal_media_metabolites:
        #        medium[key] = 0
        #    else:
        #        medium[key] = 10

        cdiff_model.enableTransportReactions()
        reactions_to_add = find_essential_reactions(cdiff_model.model)

        #reactions_to_add = cdiff_model.getFluxes().index.values
        #reactions_to_add = [cdiff_model.model.reactions.get_by_id(i) for i in reactions_to_add]
        #reactions_to_add = [reaction for reaction in reactions_to_add if "import" not in reaction.id and "Ex_" not in reaction.id and "EX_" not in reaction.id and "Trans_" not in reaction.id]


        for reaction in reactions_to_add:
            self.addCdiffReaction(cdiff_model.model, reaction.id)


        return(reactions_to_add)




    def addImportReaction(self, in_metabolite_id, reaction_id = "", reaction_name = "", ex_metabolite_id = "", ex_metabolite_name ="",
                    in_metabolite_name = "", enable_exchange = False, exchange_lim = 10):

        if reaction_id == "":
            reaction_id = "trans_" + in_metabolite_id
        if reaction_name == "":
            reaction_name = self.model.metabolites.get_by_id(in_metabolite_id).name
        if ex_metabolite_id == "":
            ex_metabolite_id = in_metabolite_id.replace("_c0", "_e0")
        if ex_metabolite_id == "":
            ex_metabolite_name = self.model.metabolites.get_by_id(in_metabolite_id).name

        if ex_metabolite_id in self.metabolite_id_list:
            ex_metabolite = self.model.metabolites.get_by_id(ex_metabolite_id)
        else:
            ex_metabolite = Metabolite(id = ex_metabolite_id, name = ex_metabolite_name)

        if in_metabolite_id in self.metabolite_id_list:
            in_metabolite = self.model.metabolites.get_by_id(in_metabolite_id)
        else:
            in_metabolite = Metabolite(id = in_metabolite_id, name = in_metabolite_name)    

        reaction = Reaction(id = reaction_id, name = reaction_name)
        reaction.add_metabolites({ex_metabolite: -1, in_metabolite: 1})
        self.model.add_reactions([reaction])

        self.model.metabolites.get_by_id(ex_metabolite_id).compartment = "e"
        self.model.metabolites.get_by_id(in_metabolite_id).compartment = "c"
        self.model = self.addCompartments()
        if enable_exchange:
            if not self.model.reactions.has_id("EX_" + ex_metabolite_id):
                self.model.add_boundary(self.model.metabolites.get_by_id(ex_metabolite_id), type = "exchange")

            medium = self.model.medium
            medium["EX_" + ex_metabolite_id] = exchange_lim
            self.model.medium = medium
        return(reaction)

    def reaction_is_reachable(self, reaction_id):
        substrates = [i.id for i in self.model.reactions.get_by_id(reaction_id).reactants]
        if len(substrates) == 0:
            return True
        if np.all([i in self.reachable_metabolite_ids for i in substrates]):
            return True
        else:
            return False
    
    def metabolite_is_reachable(self, metabolite_id):
        for reaction_id in self.reachable_reaction_ids:
            products = self.model.reactions.get_by_id(reaction_id).products
            if metabolite_id in [i.id for i in products]:
                return True
        return False

    def findUnreachables(self, reachable_reaction_ids = None):

        if reachable_reaction_ids is None:
            reachable_reaction_ids = deepcopy(self.reachable_reaction_ids)
        reachable_metabolite_ids = deepcopy(self.reachable_metabolite_ids)

        reachable_reaction_ids = set(reachable_reaction_ids)
        reachable_metabolite_ids = set(reachable_metabolite_ids)

        num_reactions_reachable = -1
        num_reactions_reachable_changed = -1
        j = 0
        #for j in range(10):
        while num_reactions_reachable != len(reachable_reaction_ids):
            num_reactions_reachable = len(reachable_reaction_ids)

            for reaction_id in reachable_reaction_ids:
                for product in self.model.reactions.get_by_id(reaction_id).products:
                    # set product to reachable
                    reachable_metabolite_ids.add(product.id)
                for reactant in self.model.reactions.get_by_id(reaction_id).reactants:
                    reachable_metabolite_ids.add(reactant.id)

            for reaction in self.model.reactions:
                if self.reaction_is_reachable(reaction.id):
                    reachable_reaction_ids.add(reaction.id)
            #j=j+1

        not_reachable_reactions = set([i.id for i in self.model.reactions]).difference(reachable_reaction_ids)
        not_reachable_metabolites = set([i.id for i in self.model.metabolites]).difference(reachable_metabolite_ids)
        #print("Reachable reactions: ", len(reachable_reaction_ids)  ) 
        #print("Reachable metabolites: ", len(reachable_metabolite_ids))
        #print("Not reachable reactions: ", len(not_reachable_reactions))
        #print("Not reachable metabolites: ", len(not_reachable_metabolites))

        return(reachable_reaction_ids, not_reachable_reactions, reachable_metabolite_ids, not_reachable_metabolites)

    def checkImportReactionPresent(self, cid):
        for reaction in self.model.reactions:
            for substrate in reaction.reactants:
                if substrate.id == cid:
                    print(substrate.id)
                    print("Import reaction present")
                    return(reaction)
        else:
            print("Import reaction not found")
            return(None)

    def checkIfCdiffReactionPresent(self, cdiff_model, reaction_id):
        reactant_ids = [self.conversions[i.id] for i in cdiff_model.model.reactions.get_by_id(reaction_id).reactants]
        product_ids = [self.conversions[i.id] for i in cdiff_model.model.reactions.get_by_id(reaction_id).products]
        for reaction in self.model.reactions:
                perc_reactant_present = np.sum([i.id in reactant_ids for i in reaction.reactants]) / len(reaction.reactants)
                perc_product_present = np.sum([i.id in product_ids for i in reaction.products]) / len(reaction.products)
                if perc_product_present == 1 and perc_product_present == 1:
                    print("Reaction already present", reaction.id, self.getReactionString(reaction.id))
                    return(True)
        print("Reaction not present: ", reaction_id)
        return(False)


    def getModelAsGraph(self):

        G = nx.DiGraph()

        # Add metabolites as nodes
        for metabolite in self.model.metabolites:
            G.add_node(metabolite.id, type='metabolite', name = metabolite.name, category = "metabolite")

        # Add reactions as nodes and connect them to metabolites
        for reaction in self.model.reactions:
            G.add_node(reaction.id, type='reaction', name = reaction.name, category = "reaction")
            if reaction.bounds[1] > 0:
                for metabolite in reaction.reactants:
                    G.add_edge(metabolite.id, reaction.id)
                for metabolite in reaction.products:
                    G.add_edge(reaction.id, metabolite.id)
            if reaction.bounds[0] < 0:
                for metabolite in reaction.reactants:
                    G.add_edge(reaction.id, metabolite.id)
                for metabolite in reaction.products:
                    G.add_edge(metabolite.id, reaction.id)    

        self.G = G
        return(G)

    def getModelAsReactionGraph(self):
        G = self.getModelAsGraph()
        # Initialize a new directed graph for reactions
        reaction_graph = nx.DiGraph()

        # Separate the reactions and compounds from the original graph
        reactions = [n for n, attr in G.nodes(data=True) if attr['type'] == 'reaction']
        compounds = [n for n, attr in G.nodes(data=True) if attr['type'] == 'metabolite']

        # Add all reactions as nodes in the new graph
        reaction_graph.add_nodes_from((n, G.nodes[n]) for n in reactions)

        # Loop through each compound and find the reactions connected to it
        for compound in compounds:
            # Find all reactions connected to this compound
            predecessors = list(G.predecessors(compound))  # Incoming edges (reactions -> compound)
            successors = list(G.successors(compound))      # Outgoing edges (compound -> reactions)

            # Add edges between reactions if they share the same compound
            for pred in predecessors:
                for succ in successors:
                    # Add a directed edge from the predecessor reaction to the successor reaction
                    reaction_graph.add_edge(pred, succ, title = self.model.metabolites.get_by_id(compound).name)

        return(reaction_graph)

    def calculateMinimalMedia(self, percent_max_growth = 1, minimize_components = False):
        # 1. Load gapfilled model w/ added transport reactions from C. diff


        # 2. Convert all bi-directional reactions to be single uni-directional reactions:
        #model = mm.convertReactionsUnidirectional(model)

        # 4. Calculate minimal media formulation using cobra

        self.model = self.addCompartments()
        max_growth = percent_max_growth * self.model.slim_optimize()
        minimal_media_metabolites = minimal_medium(self.model, max_growth, minimize_components = minimize_components).index.values
        #minimal_media_metabolites = [i.replace("EX_", "") for i in minimal_media_metabolites]
        print("Number of metabolites in minimal media: ", len(minimal_media_metabolites))
        
        return(minimal_media_metabolites)

    def setMinimalMedia(self, minimal_media_metabolites):
        medium = self.model.medium
        for met_id in self.model.medium.keys():
            medium[met_id] = 0
        for met_id in minimal_media_metabolites:
            media_id = "EX_" + met_id.replace("_rev", "")
            if media_id in [i.id for i in self.model.boundary]:
                medium[media_id] = 10
        self.model.medium = medium
        return(minimal_media_metabolites)

    def setCompleteMedia(self):
        medium = self.model.medium
        for met_id in self.model.medium.keys():
            medium[met_id] = 10
        self.model.medium = medium

    def convertReactionsUnidirectional(self):
        new_reactions = []

        # Iterate through each reaction in the model
        for reaction in self.model.reactions:
            if reaction.lower_bound < 0:
                # Create a new reaction
                new_reaction = Reaction(reaction.id + "_rev", name = reaction.name + "_rev")  # Name the new reaction
                new_reaction.lower_bound = 0  # Set appropriate bounds (you might want to adjust this)
                new_reaction.upper_bound = -reaction.lower_bound  # You may choose to set it to a similar upper bound
                new_reaction.gene_reaction_rule = reaction.gene_reaction_rule
                
                # Switch the reactants and products
                new_reaction.add_metabolites({metabolite: -coef for metabolite, coef in zip(reaction.reactants, [reaction.metabolites[i] for i in reaction.reactants])})
                new_reaction.add_metabolites({metabolite: -coef for metabolite, coef in zip(reaction.products, [reaction.metabolites[i] for i in reaction.products])})

                # Add the new reaction to the model
                self.model.reactions.get_by_id(reaction.id).lower_bound = 0
                new_reactions.append(new_reaction)
        self.model.add_reactions(new_reactions)

    def getReactionString(self, reaction_id):
        reaction = self.model.reactions.get_by_id(reaction_id)
        reactant_names = [i.name.replace(" [c0]", "") for i in reaction.reactants]
        product_names = [i.name.replace(" [c0]", "") for i in reaction.products]
        reaction_str = ' + '.join(reactant_names) + "=>" + " + ".join(product_names)
        return(reaction_str)

    def saveModel(self, filepath):
        cobra.io.save_json_model(self.model, filepath)

    def plotReactionNetwork(self, G, reaction_ids, reaction_to_color = None, edge_derivatives = []):
        net = Network(notebook=True, directed = True, height="1000px", width="1600px")
        subgraph = G.subgraph(list(reaction_ids))
        # Add nodes from networkx to pyvis, preserving 'name' attribute as label
        for node, data in subgraph.nodes(data=True):
            if reaction_to_color is not None:
                net.add_node(node, label=data['name'],
                     color = "pink",
                    title = self.getReactionString(node),
                    size = 14)  # Use 'name' attribute as label
            net.add_node(node, label=data['name'],
                        title = self.getReactionString(node),
                        size = 14)  # Use 'name' attribute as label
            
        for source, target, data in subgraph.edges(data = True):
            compound_name = data['title']
            edge_width = 0.6
            color_list = ["rgba(18, 188, 226, 0.5)",  "rgba(237, 106, 90, 0.5)", "rgba(236, 203, 70, 0.5)"]
            color="rgba(0,0,0,0.2)"
            edge_width = 1
            for i, derivative_list in enumerate(edge_derivatives):
                if np.sum([i in compound_name for i in derivative_list]) > 0:
                    color = color_list[i]
                    edge_width = 5
            net.add_edge(source, target, color=color, title = compound_name, width = edge_width)
        return(net)


    def getModelDiffs(self, model_orig):
        reactions_added = [i for i in self.model.reactions if i not in model_orig.model.reactions]
        reactions_added = [i for i in reactions_added if "SK_" not in i.id and "test" not in i.id]
        reactions_added = [i for i in reactions_added if "EX_" not in i.id]
        ids = [i.id for i in reactions_added]
        names = [i.name for i in reactions_added]
        reaction_strings = [self.getReactionString(i.id) for i in reactions_added]
        gene_rules  = [self.model.reactions.get_by_id(i.id).gene_reaction_rule for i in reactions_added]
        df = pd.DataFrame({'id': ids, 'name': names, 'reaction_string': reaction_strings, 'gene_rule': gene_rules})
        df.index = df['id']
        return(df)



    def growth_function_binary(self, selected_metabolites_bool):

        # 1. set model medium to minimal medium
        self.model.medium = deepcopy(self.minimal_medium)
        #print("minimal media is: ", len(self.model.medium))

        # 2. fetch possible metabolite search space
        medium_components = self.medium_search_space

        # 3. get names of the metabolites selected by the optimizer
        selected_met_ids = np.array(medium_components)[[i == 1 for i in selected_metabolites_bool]]
        #print(selected_metabolites_bool)
        #print(selected_met_ids)

        # 4. Increase each by a big step
        medium = self.model.medium
        for met_id in selected_met_ids:
            medium[met_id] = 1000
        self.model.medium = medium

        # 5. return the growth value under those conditions
        growth = self.model.optimize().objective_value

        return -growth  # We use a negative because gp_minimize minimizes the objective


    def runOptimize(self, medium_components, n_calls = 25, seed = 0, verbose = False):
        self.medium_search_space = medium_components
        n_metabolites = len(self.medium_search_space)
        
        search_space = [Categorical([0, 1], name=self.medium_search_space[i], prior = [0.90, 0.1]) for i in range(n_metabolites)]
        
        result = gp_minimize(
            func=self.growth_function_binary,              # The objective function to minimize
            dimensions=search_space,     # The search space
            n_calls=n_calls,                  # Number of function evaluations
            random_state=seed,               # For reproducibility
            verbose = verbose
        )

        best_selection = result.x
        best_selection = [i == 1 for i in best_selection]
        selected_medium_components = np.array(medium_components)[best_selection]

        growth_rate = -self.growth_function(best_selection)
        return(selected_medium_components, best_selection, growth_rate, result)


    def runOptimizeLinear(self, medium_components, n_calls = 25, seed = 0, verbose = False):

        self.medium_search_space = medium_components
        n_metabolites = len(self.medium_search_space)
        
        search_space = [Categorical([0, 1], name=self.medium_search_space[i], prior = [0.90, 0.1]) for i in range(n_metabolites)]
        
        result = minimize(
            fun=self.growth_function,              # The objective function to minimize
            dimensions=search_space,     # The search space
            n_calls=n_calls,                  # Number of function evaluations
            random_state=seed,               # For reproducibility
            verbose = verbose
        )

        best_selection = result.x
        best_selection = [i == 1 for i in best_selection]
        selected_medium_components = np.array(medium_components)[best_selection]

        growth_rate = -self.growth_function(best_selection)
        return(selected_medium_components, best_selection, growth_rate, result)


    def growth_function(self, metabolite_concentrations):

        # 1. set model medium to minimal medium
        self.model.medium = deepcopy(self.minimal_medium)
        #print("minimal media is: ", len(self.model.medium))

        # 2. fetch possible metabolite search space
        medium_components = self.medium_search_space

        # 3. get names of the metabolites selected by the optimizer
        #selected_met_ids = np.array(medium_components)[[i == 1 for i in selected_metabolites_bool]]

        # 4. Increase each by a big step
        medium = self.model.medium
        for i, met_id in enumerate(medium_components):
            medium[met_id] = metabolite_concentrations[i]
        self.model.medium = medium

        # 5. return the growth value under those conditions
        growth = self.model.optimize().objective_value
        
        reg_strength = 0.4
        regularize = np.sum(np.abs(metabolite_concentrations))
        return -growth + (reg_strength * regularize)  # We use a negative because minimize minimizes the objective


    def getReactionCorrelations(self, selected_met_ids):
        self.model.medium = deepcopy(self.minimal_medium)
        fluxes_list = []
        concentrations = [1, 2]
        for conc in concentrations:
            medium = deepcopy(self.minimal_medium)
            for met_id in selected_met_ids:
                medium["EX_" + met_id] = conc
                model.model.medium = medium
            print(len(model.model.medium))
            fluxes = model.model.optimize().fluxes
            print( model.model.optimize().objective_value)
            fluxes_list.append(fluxes)
            model.model.medium = deepcopy(self.minimal_medium)

        fluxes_df = pd.concat(fluxes_list, axis=1)

        reaction_ids_correlated = []
        stats = []
        pvals = []
        for i in range(fluxes_df.shape[0]):
            stat = pearsonr(fluxes_df.iloc[i, :].values, concentrations)[0]
            pval = pearsonr(fluxes_df.iloc[i, :].values, concentrations)[1]
            if pval < 1:
                reaction_ids_correlated.append(fluxes_df.index.values[i])
                stats.append(stat)
                pvals.append(pval)
        stat_df = pd.DataFrame({'stats': stats, 'pvals': pvals}, index = reaction_ids_correlated)

        return(stat_df)

    def value_to_color(value):
        # Normalize the value to be between 0 and 1 for colormap
        normalized_value = (value + 1) / 2

        # Use a colormap (RdBu is a good choice for blue to red)
        cmap = plt.get_cmap('RdBu')

        # Get the RGBA color from the colormap
        rgba_color = cmap(normalized_value)

        # Convert RGBA to hex format
        rgba_string = 'rgba({:.0f},{:.0f},{:.0f},{:.2f})'.format(
            rgba_color[0] * 255,  # Red
            rgba_color[1] * 255,  # Green
            rgba_color[2] * 255,  # Blue
            0.2                  # Alpha (transparency)
        )
        
        return rgba_string

    def colorGraph(graph):
        # Add nodes from networkx to pyvis, preserving 'name' attribute as label
        for node, data in graph.nodes(data=True):
            if 'correlation' not in data:
                data['correlation'] = 1
            color = model.value_to_color(data['correlation'])
            if node == 'rxn01517_c0_rev':
                color = "pink"
            net.add_node(node, label=data['name'],
                         color = color,
                        title = model.getReactionString(node),
                        size = 14)  # Use 'name' attribute as label
            
            
        # Add edges from networkx to pyvis
        for source, target, data in subgraph.edges(data = True):
            compound_name = data['title']
            edge_width = 1
            if np.sum([i in compound_name for i in nad_derivative]) > 0:
                color = "rgba(237, 106, 90, 0.5)"
                edge_width = 5
            elif np.sum([i in compound_name for i in atp_derivative]) > 0:
                color = "rgba(236, 203, 70, 0.5)"
                edge_width = 5
            else:
                color="rgba(0,0,0,0.3)"
            net.add_edge(source, target, color=color, title = compound_name, width = edge_width)



class database:
    def __init__(self):
        self.db = pd.DataFrame()
        
    def addInfo(self, data, col_name):
        #data should be a rxn: value dictionary
        self.db[col_name] = self.db.index.map(data)
        self.db[col_name] = self.db[col_name].replace({np.nan: None})
        return(self.db)
    
    #def addInfoNewKey(self, key_col_name, data, col_name):
        
    def getModelSeedReactionDB(self):
        # Read in both sources (both from github) of name to reaction data
        modelseed_reactions = pd.read_csv("../data/model_seed_database/modelSEED_reactions.tsv", sep = "\t")
        modelseed_reactions
        seed_db = pd.read_csv("../data/model_seed_database/ModelSEED_Subsystems.tsv", sep = "\t")
        seed_db.head()
        seed_db.columns = ["Class", "Sub-class", "Cateogry", "name", "id"]

        rxn_to_name = pd.concat([modelseed_reactions, seed_db])
        rxn_to_name = rxn_to_name.set_index('id')

        return(rxn_to_name)
    
    def getRxnForEC(self, ec):
        rxn_ids = self.db[self.db['EC_from_db'] == ec].index.values
        return(rxn_ids)

    def getNameForEC(self, ec):
        names = self.db[self.db['EC_from_db'] == ec].name.values
        return(names)

    def addTranscriptomicData(self, ec_numbers):
        # add in transcriptomics data
        tmp = self.db
        tmp['rxn_id'] = tmp.index.values
        tmp = tmp.set_index("EC_from_db")
        tmp['present_in_pbi_transcript'] = ["" for i in range(tmp.shape[0])]

        for ec in ec_numbers:
            tmp.loc[ec, "present_in_pbi_transcript"] = ec
        tmp['EC_from_db'] = tmp.index

        self.db = tmp.set_index('rxn_id')
        self.db = self.db[self.db.index.notnull()]
        return(self)

def createDB(import_model, db):
    ecs_from_model = extract_ec_numbers_from_json(import_model)
    db.addInfo(ecs_from_model, "EC_from_modelseed").head(2)

    # read in rxn to ec from http://localhost:8888/notebooks/aidan/sandbox/compare_model_w_imported_model.ipynb

    df = pd.read_csv("../data/model_seed_database/rxn_to_ec.csv", index_col = 0)
    df['EC'] = [i.replace("EC-", "") for i in df['EC']]
    db.addInfo(df['EC'].to_dict(), col_name = "EC_from_db").head(2)

    names_from_model = extract_names_from_model(import_model)
    db.addInfo(names_from_model, "name_from_modelseed").head(2)

    # add all reaction aliases from database
    df = pd.read_csv("../data/model_seed_database/Unique_ModelSEED_Reaction_Names.txt", sep = "\t")
    df['External ID'] = df['External ID'].astype(str).replace('None', '')

    df = df.groupby("ModelSEED ID")['External ID'].apply(lambda x: ', '.join(sorted(set(x)))).reset_index().set_index("ModelSEED ID")
    db.addInfo(df["External ID"].to_dict(), "names_from_seeddb").head(2)

    gene_rules_from_model = extract_generule_from_model(import_model)
    db.addInfo(gene_rules_from_model, 'pbi_gene_rule').head(3)
    return(db)



# Functions to read json files from KBASE:
def extract_ec_numbers_from_reaction(reaction):
    ec_numbers = []
    # Traverse the 'modelReactionProteins' list
    for protein in reaction.get('modelReactionProteins', []):
        for subunit in protein.get('modelReactionProteinSubunits', []):
            note = subunit.get('note', '')
            # Use regular expression to find EC numbers
            matches = re.findall(r'EC (\d+\.\d+\.\d+\.\d+)', note)
            ec_numbers.extend(matches)
    return ec_numbers

def extract_ec_numbers_from_json(import_model):
    ecs_from_model = {}
    for reaction in import_model.get('modelreactions'):
        reaction_id = reaction.get('id').replace("_c0", "").replace("_e0", "").replace("_z0", "")
        ecs_from_model[reaction_id] = extract_ec_numbers_from_reaction(reaction)
    return(ecs_from_model)

def extract_names_from_reaction(reaction):
    names = []
    for protein in reaction.get('modelReactionProteins', []):
        for subunit in protein.get('modelReactionProteinSubunits', []):
            note = subunit.get('note')
            names.append(note)
    return names

def extract_names_from_model(import_model):
    names_from_model = {}
    for reaction in import_model.get('modelreactions'):
        reaction_id = reaction.get('id').replace("_c0", "").replace("_e0", "").replace("_z0", "")
        names_from_model[reaction_id] = extract_names_from_reaction(reaction)
    return(names_from_model)

# get CDS gene rule from model
def extract_generule_from_reaction(reaction):
    patterns = []
    for protein in reaction.get('modelReactionProteins', []):
        for subunit in protein.get('modelReactionProteinSubunits', []):
            feature_refs = subunit.get('feature_refs', [])
            # Extract pattern using regular expressions
            for ref in feature_refs:
                matches = re.findall(r'CDS\.\d+', ref)
                patterns.extend(matches)
    return patterns

def extract_generule_from_model(import_model):

    gene_rules_from_model = {}
    for reaction in import_model.get('modelreactions'):
        reaction_id = reaction.get('id').replace("_c0", "").replace("_e0", "").replace("_z0", "")
        gene_rules_from_model[reaction_id] = list(np.unique(extract_generule_from_reaction(reaction)))
    return(gene_rules_from_model)

def cleanECs(ecs):
    ecs = [i for i in ecs if '-' not in i and all(part.isdigit() for part in i.split('.'))]
    return(ecs)

def buildStoichMatrix(import_model):

    # make a reaction matrix from the imported model
    reaction_ids = []
    for reaction in import_model['modelreactions']:
        reaction_ids.append(reaction['id'])

    compounds_list = []
    for compound in import_model['modelcompounds']:
        compounds_list.append(compound['id'])

    s_import = pd.DataFrame(index = reaction_ids, columns = compounds_list)

    s_import
    for reaction in import_model['modelreactions']:
        reaction_id = reaction['id']
        compounds = [i['modelcompound_ref'].split("/")[3] for i in reaction['modelReactionReagents']]
        for i, compound in enumerate(compounds):
            s_import.loc[reaction_id, compound] = reaction['modelReactionReagents'][i]['coefficient']

    s_curated = deepcopy(s_import)
    s_curated = s_curated.fillna(0)
    #s_curated.head()

    s_curated.columns = [i.split("_")[0] for i in s_curated.columns]
    print(s_curated.shape)
    
    return(s_curated)

def getInfoFromKBASE(filename):
    f = open(filename)
    import_model = json.load(f)
    print(len(import_model['modelreactions']))
    f.close()

    # 2. Read in conversions to bigg
    conversion = pd.read_csv("../data/model_seed_database/bigg_to_seed.tsv", sep = "\t", index_col = 0)
    seed_to_bigg = {}
    bigg_to_seed = {}
    for i in range(conversion.shape[0]):
        bigg = conversion.index.values[i]
        seed = conversion.seed.values[i]
        seed_to_bigg[seed] = bigg
        bigg_to_seed[bigg] = seed

    # 3. Make a matrix S
    reaction_ids = []
    reaction_names = {}
    lower_bounds = {}
    upper_bounds = {}
    for reaction in import_model['modelreactions']:
        reaction_ids.append(reaction['id'])
        reaction_names[reaction['id']] = reaction['name']
        lower_bounds[reaction['id']] = -reaction['maxrevflux']
        upper_bounds[reaction['id']] = reaction['maxforflux']

    compounds_list = []
    compound_names = {}
    for compound in import_model['modelcompounds']:
        compounds_list.append(compound['id'])
        compound_names[compound['id']] = compound['name']

    s_import = pd.DataFrame(index = reaction_ids, columns = compounds_list)

    for reaction in import_model['modelreactions']:
        reaction_id = reaction['id']
        compounds = [i['modelcompound_ref'].split("/")[3] for i in reaction['modelReactionReagents']]
        for i, compound in enumerate(compounds):
            s_import.loc[reaction_id, compound] = reaction['modelReactionReagents'][i]['coefficient']


    s_auto = deepcopy(s_import)
    s_auto = s_auto.fillna(0)
    print(s_auto.shape)
    s_auto.head()
    return(s_auto, reaction_names, compound_names, lower_bounds, upper_bounds)

def getCobraModel(db, s, reaction_names, compound_names, lower_bound_dict, upper_bound_dict):
    model = Model('example_model')
    
    #1.  define some metabolites
    metabolite_dict = {}
    for compound_id in s.columns.values:
        metabolite_dict[compound_id] = Metabolite(id = compound_id, name = compound_names[compound_id]) # add old name as the formula
    
    for reaction_id in s.index.values:
        reaction = Reaction(id = reaction_id,
                name = reaction_names[reaction_id],
                lower_bound = lower_bound_dict[reaction_id],
                upper_bound = upper_bound_dict[reaction_id])

        reaction_compounds = s.loc[reaction_id][s.loc[reaction_id] != 0]
        compound_dict = {}
        for c_id in reaction_compounds.index.values:
            if np.sum([c_id == i for i in reaction_compounds.index.values]) > 1:
                compound_dict[metabolite_dict[c_id]] = reaction_compounds[c_id][0]
                compound_dict[metabolite_dict[c_id]] = reaction_compounds[c_id][1]
            else:
                compound_dict[metabolite_dict[c_id]] = reaction_compounds[c_id]
        reaction.add_metabolites(compound_dict)
        model.add_reactions([reaction])
    return(model)

#########################

def addNewRXNs(db, rxn_to_add, model):
#TODO: address reversibility
#TODO: address issue with _c vs. _e if we're adding the reaction
    rxns_in_model = [i.id for i in model.reactions]
    rxns_in_model = delete_component(rxns_in_model)

    new_metabolites = {}

    for rxn_id in rxn_to_add:
        if rxn_id not in rxns_in_model:
            rxn_id_new = rxn_id + "_c0"
            if not model.reactions.has_id(rxn_id_new):
                print("Attempting to add: ", rxn_id_new)
                compound_dict = {}
                num_entries = sum([i == rxn_id for i in db.db.index.values])
                reaction_db = None
                if  num_entries > 1:
                    reaction_db =  db.db.loc[rxn_id, :]
                    reaction_db = reaction_db.iloc[0, :]
                    stoich_str = reaction_db['stoichiometry']
                    #stoich_str = db.db.loc[rxn_transcript_unique[2], "stoichiometry"].values[0]
                else:
                    reaction_db =  db.db.loc[rxn_id, :]
                    stoich_str = reaction_db['stoichiometry']
                    #stoich_str = db.db.loc[rxn_id, "stoichiometry"]

                items = stoich_str.split(";")
                try:
                    for item in items:
                        coefficient = item.split(":")[0]
                        compound_id = item.split(":")[1] + "_c0" 
                        compound_name = item.split(":")[4]
                        if not model.metabolites.has_id(compound_id):
                            new_metabolites[compound_id] = Metabolite(id = compound_id, name = compound_name)
                            model.add_metabolites(new_metabolites[compound_id])
                        compound_dict[compound_id] = int(coefficient)

                    # make reaction

                    reaction_name = str(reaction_db['name'])

                    reaction = Reaction(id = rxn_id_new,
                         name = reaction_name,
                         lower_bound = -100,
                         upper_bound = 100)

                    if not model.reactions.has_id(rxn_id_new):
                        model.add_reactions([reaction])

                        reaction.add_metabolites(compound_dict)
                except:
                    print("Error adding: ", rxn_id_new)
        else:
            print(rxn_id, " already present")
    return(model)




def delete_component(reactions):
    return([i.replace("_c0", "").replace("_z0", "").replace("_e0", "") for i in reactions])

def delete_sul(reactions):
    return([i for i in reactions if i != 'sul00010'])






def addBiomassReaction(model, import_model, biomass_reaction_id):
    reaction = Reaction(id = biomass_reaction_id, name = "biomass", lower_bound = 0, upper_bound = 1000)

    metabolite_dict = {}
    for compound in import_model['biomasses'][0]['biomasscompounds']:
        coefficient = compound['coefficient']
        cid = compound['modelcompound_ref'].split("/")[3]
        if not model.metabolites.has_id(cid):
            print(cid, " was missing, now added")
            met = Metabolite(id = cid, name = "")
        else:
            met = model.metabolites.get_by_id(cid)
        metabolite_dict[met] = coefficient
    reaction.add_metabolites(metabolite_dict)
    model.add_reactions([reaction])
    return(model)




# add an atp sink reaction
def addATPSinkObjective(model):
    atp = "cpd00002_c0"

    reaction = Reaction(id = "ATP_sink",
            name = "ATP objective",
            lower_bound = 0,
            upper_bound = 1000)
    reaction.add_metabolites({model.metabolites.get_by_id(atp): -1})
    model.add_reactions([reaction])
    model.reactions.get_by_id("ATP_sink")

    objective_list = ["ATP_sink"]
    objective_dict = {}
    for objective in objective_list:
        objective_reaction = model.reactions.get_by_id(objective)
        objective_dict[objective_reaction] = 1 / len(objective_list) # assumes all equally important
    model.objective = objective_dict
    model.reactions.get_by_id(objective).upper_bound=1000
    return(model)



def addDummyDatabaseEntry(id_, db):
    atp_reaction_entry = [None for i in range(db.shape[1])]
    atp_reaction_entry = pd.DataFrame([None for i in range(db.shape[1])], index = db.columns.values).T
    atp_reaction_entry.index = [id_]
    db = pd.concat([db, atp_reaction_entry])
    return(db)


def addMediaBoundaries(model, db):

    media = pd.read_csv("/home/christine/Partners HealthCare Dropbox/Christine Tataru/metabolic_modeling/data/nat_chem_bio_paper/c_scindens_mm.tsv", sep = "\t", index_col = 0)
    media.index = [i + "_e0" for i in media.index.values]
    model.metabolites.has_id
    media_dict = {}
    for compound in media.index.values:
        print(compound)
        if model.metabolites.has_id(compound):
            new_id = 'EX_' + compound
            if not model.reactions.has_id(new_id):
                model.add_boundary(model.metabolites.get_by_id(compound), type = "exchange")
                media_dict[new_id] = media.loc[compound, "concentration"]
                db = addDummyDatabaseEntry(new_id.replace("_e0", ""), db)

    model.medium = media_dict
    return(model, db)




def findManualCompoundConversions(cdiff_model, reaction_id):
    db_compound = pd.read_csv("../data/model_seed_database/compounds.tsv", sep = "\t", index_col = 0)

    sp = cdiff_model.reactions.get_by_id(reaction_id)
    ids = [i.id for i in sp.metabolites]
    names = [cdiff_model.metabolites.get_by_id(i).name for i in ids]
    for name in names:
        print(name)
        keep = [name in str(alias) for alias in db_compound['aliases'].values]
        display(db_compound.loc[keep, :].iloc[0:4, :])










def testModelBiomass(model, objective = "bio_kbase"):
    model = deepcopy(model)
    model = convertReactionsUnidirectional(model)
    
    # make all transport reactions accessible
    transport_reactions = getTransportReactions(model)
    model = addCompartments(model)
    model, boundary_to_transport = enableTransportReactions(model, transport_reactions)

    #optim = model.optimize()
    optim = flux_variability_analysis(model, fraction_of_optimum=0.95)
    return(optim, boundary_to_transport)





def getMissingComponents(model, reaction_id, reachable_reaction_ids):
    print(reaction_id)
    missing_components = []
    reachable_df = []
    reactants = []
    for reactant in model.reactions.get_by_id(reaction_id).reactants:
        reachable = metabolite_is_reachable(model, reactant.id, reachable_reaction_ids)
        reachable_df.append(reachable)
        reactants.append(reactant.name)
        #print(reactant.name, ": ", reachable)
        if not reachable:
            missing_components.append(reactant.id)
    reachable_df = pd.DataFrame(reachable_df, index = reactants)
    return(missing_components, reachable_df)


def addTranscriptomicReactions(model, db):
    ecs_to_add = pd.read_csv("../data/nat_chem_bio_paper/transcriptomics/rast/ecs_unique.txt", header = None).iloc[:, 0].values
    db_ec = db.db.set_index('EC_from_db')
    db_ec['rxn'] = db.db.index.values

    rxn_ids = db_ec.loc[[i for i in ecs_to_add if i in db_ec.index.values], "rxn"].values
    rxn_names = db_ec.loc[[i for i in ecs_to_add if i in db_ec.index.values], "name"].values
    rxn_names = [str(i) for i in rxn_names]
    rxn_names = np.unique(rxn_names)
    model = addNewRXNs(db, rxn_to_add = rxn_ids, model = model)
    return(model, db)



def get_db_compound(compounds_filepath):
    db_compound = pd.read_csv(compounds_filepath, sep = "\t", index_col = 0)
    db_compound.index = [i + "_c0" for i in db_compound.index.values]
    db_compound_copy = deepcopy(db_compound)
    db_compound_copy.index = [i.replace("_c0", "_e0") for i in db_compound.index.values]
    db_compound = pd.concat([db_compound, db_compound_copy])
    return(db_compound)

def getConversionFromCdiff(compounds_filepath, curated_file_path, cdiff_model_filepath):
    curated_to_seed_compounds = pd.read_csv(curated_file_path)
    curated_to_seed_compounds

    # conversion method 1
    conversions = curated_to_seed_compounds.groupby('orig_IDs')['new_IDs'].agg(''.join).to_dict()
    for key in conversions.keys():
        if "_c" in key:
            conversions[key] = conversions[key] + "_c0"
        if "_e" in key:
            conversions[key] = conversions[key] + "_e0"

    print("Length of conversions after 1st method: ", len(conversions))
    # conversion method 2
    compounds = pd.read_csv(compounds_filepath, sep = "\t").set_index("name")
    compounds.index = [str(i).lower() for i in compounds.index.values]
    compounds
    f = open(cdiff_model_filepath)
    json_model = json.load(f)
    f.close()
    cdiff_metabolites = {}
    for met in json_model['metabolites']:
        cdiff_metabolites[met['id']] = met['name']
        
    for key, name in cdiff_metabolites.items():
        name = name.lower()
        if not name in compounds.index.values:
            name = name.replace("2+", "")
        if name in compounds.index.values:
            ending = ""
            if "_c" in key:
                ending = "_c0"
            if "_e" in key:
                ending = "_e0"
            if key not in conversions:
                conversions[key] = compounds.loc[name, "id"] + ending
    print("Length of conversions after 2nd method: ", len(conversions))

    #conversion method 3
    conversions.update({'h2o': 'cpd00001',
                       'pi': 'cpd00009',
                       'tagd': 'cpd00437',
                       'ppa' : 'cpd00141',
                       'ppat' : 'cpd04099',
                       'gluL': 'cpd00023',
                       'asnL': 'cpd00132',
                       'glyb': 'cpd00540', 
                       'aspL': 'cpd00041',
                        'ptrc' : 'cpd00118',
                        'metD' : 'cpd00637',
                        'chol' : 'cpd00098',
                        'fe2' : 'cpd00021',
                        'c1ala' : "cpd12878",
                        'raf' : "cpd00382", # double check with lynn
                        'udpagaepyr': 'cpd02820', # check that this should be a transport reaction - big molecule
                       'hxan': 'cpd00226',
                        'hco3': 'cpd00242',
                        '4abz': 'cpd00443',
                        'pntoR': 'cpd00644',
                        'nh3': 'cpd00013',
                        'bio': 'cpd00104',
                        '2o6pamg2d': 'cpd21036', # double check; only option already present in model
                       'sor': 'cpd00588',
                        'mevR': 'cpd00332', #double check R vs. S
                        'isobuta' : 'cpd01711', #double check
                        'ival': 'cpd05178', #double check
                        '4s45dhp23do': "cpd08638",
                        'fum' : 'cpd00106',
                        '23diap' : 'cpd03828',
                        'pyin' : 'cpd00263',
                        'pval' : 'cpd00263',
                        '2aepat' : 'cpd02233',
                        'ibtol' : 'cpd10408',
                        '2mbut' : 'cpd19585',
                        'nac'  : 'cpd00218',
                        'etoa' : 'cpd00162',
                        '4mpo' : 'cpd01042',
                        'tmam' : 'cpd00441',
                        'glcn' : 'cpd00222',
                        'celb' : 'cpd00158',
                        'lacS' : 'cpd00159',
                        'malS' : 'cpd00130',
                        '2dh3dg' : 'cpd00176', #double check
                        'xan' : 'cpd00309',
                        'h2s' : 'cpd00239',
                        'gtol': 'cpd01171',
                        'tre': 'cpd00794',
                        'mnl' : 'cpd00314',
                        'dachi' : 'cpd01157',
                        'tgt' : 'cpd00589',
                        'fuc' : 'cpd01186',
                        'ca2' : 'cpd00063',
                        'acoa' : 'cpd00022',
                        'citrL' : 'cpd00274',
                        'utp': 'cpd00062',
                        'atp' : 'cpd00002',
                        'adp' : 'cpd00008',
                        'nad' : 'cpd00003',
                        'nadh' : 'cpd00004',
                        'nadph' : 'cpd00005',
                        'nadp' : 'cpd00006',
                        'ctp' : 'cpd00052',
                        'gtp' : 'cpd00038',
                        'g3p' : 'cpd00080',
                        'nac': 'cpd00218',
                        'nacmnc' : 'cpd28952',
                        'nad' : 'cpd00003',
                        'hco3' : 'cpd00242',
                        'gal' : 'cpd00709',
                        'palmacd' : 'cpd00214',
                        'myrstcoa' : 'cpd01695',
                        'serL' : 'cpd00054',
                        'gly' : 'cpd00033',
                        'ileL' : 'cpd00322',
                        'duri' : 'cpd00412',
                        'udpamr' : 'cpd01757',
                        'aga_c': 'cpd00122',
                        'dgsn' : 'cpd00277',
                        'din': 'cpd03279',
                        'cytd': 'cpd00367',
                        'proL' : 'cpd00129',
                        'cysL': 'cpd00084',
                        'leuL' : 'cpd00107',
                        'valL' : 'cpd00156',
                        'ileL' : 'cpd00322',
                        'thrL' : 'cpd00161',
                        'gly': 'cpd00033',
                        'metL': 'cpd00060',
                        'argL': 'cpd00051',
                        'alaD': 'cpd00117'

                       })
    print("Length of conversions after 3rd method: ", len(conversions))
    conversions_copy = deepcopy(conversions)
    for key in conversions_copy.keys():
        if "_c" not in key and "_e" not in key:
            conversions[key + "_c"] = conversions[key] + "_c0"
            conversions[key + "_e"] = conversions[key] + "_e0"
            conversions.pop(key)
    return(conversions)






def visualize_n_hops(G, start_node, n, reachable_reaction_ids, reachable_metabolite_ids):
    conversions = getConversionFromCdiff()
    met_ids = ['atp_c', 'adp_c', 'nad_c', 'nadh_c', 'nadph_c', 'nadp_c', 'ppi_c', 'h2o_c', 'pi_c', 'utp_c', 'ctp_c',
              'gtp_c', 'coa_c', 'fad_c', 'acoa_c']
    met_ids = [conversions[i] for i in met_ids]


    # Get nodes within n hops using BFS
    nodes_within_n_hops = set()
    queue = [(start_node, 0)]  # (current node, current hop level)
    
    while queue:
        current_node, current_hop = queue.pop(0)
        if current_hop > n:
            continue
        nodes_within_n_hops.add(current_node)
        
        # Get all outgoing neighbors
        outgoing_neighbors = set([v for u, v in G.out_edges(current_node)])
        # Get all incoming neighbors
        incoming_neighbors = set([u for u, v in G.in_edges(current_node)])
        # Combine both to get all neighbors
        all_neighbors = outgoing_neighbors.union(incoming_neighbors)

        print(current_node)
        for neighbor in all_neighbors:
            if neighbor not in nodes_within_n_hops:
                if neighbor not in reachable_reaction_ids or neighbor not in reachable_metabolite_ids:
                    if neighbor not in met_ids:
                        if neighbor not in ['Bio_lip', 'rxn05296_c0']:
                            queue.append((neighbor, current_hop + 1))

    # Include both incoming and outgoing edges for the nodes within n hops
    for node in list(nodes_within_n_hops):
        # Add outgoing edges
        for neighbor in G.neighbors(node):
            nodes_within_n_hops.add(neighbor)
        # Add incoming edges
        for predecessor in G.predecessors(node):
            nodes_within_n_hops.add(predecessor)

    # Create a subgraph containing only the relevant nodes and their edges
    subgraph = G.subgraph(nodes_within_n_hops)

    # Define a color map based on the 'category' attribute
    categories = {node: G.nodes[node]['category'] for node in subgraph.nodes()}
    unique_categories = list(set(categories.values()))
    color_map = {category: plt.cm.tab10(i) for i, category in enumerate(unique_categories)}
    color_map2 = {'True': "blue", 'False':"red"}

    # Assign colors to nodes
    #node_colors = [color_map[categories[node]] for node in subgraph.nodes()]
    reachability = [i in reachable_metabolite_ids or i in reachable_reaction_ids for i in subgraph.nodes()]
    node_colors = [color_map2[str(i)] for i in reachability]

    # Draw the graph with increased spacing between nodes
    #pos = nx.spring_layout(subgraph, k=1.5)  # Increase k for more spacing


    pos = nx.nx_agraph.graphviz_layout(subgraph, prog="dot")


    nx.draw(subgraph, pos, with_labels=True, node_color=node_colors, node_size=2000, font_size=10, arrows=True)
    plt.title(f'Graph from {start_node} with {n} hops (colored by category)')
    plt.show()



def addAllCdiffBiomassReactions(model, cdiff_model):

    # Used to make the conversion dictionaries originally
    # #mm.findManualCompoundConversions(cdiff_model, "Bio_CLP")

    reaction_id = "Bio_SPs"
    db_compound = pd.read_csv("../data/model_seed_database/compounds.tsv", sep = "\t", index_col = 0)
    db_compound.index = [i + "_c0" for i in db_compound.index.values]
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_sp)

    reaction_id = "Bio_lip"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_lip)

    reaction_id = "Bio_prot"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_prot)

    reaction_id = "Bio_CLP"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_clp)

    reaction_id = "Bio_DNA"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_dna)

    reaction_id = "Bio_RNA"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_rna)

    reaction_id = "Bio_CW"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio_cw)

    reaction_id = "bio1"
    model, db_compound, reaction = addCdiffReaction(model, cdiff_model, reaction_id, db_compound, conversions_bio1)

    return(model, db_compound)


conversions_transport = {
    'adp_c': "cpd00008_c0",
    'atp_c': "cpd00002_c0",
    'h2o_c': 'cpd00001_c0',
    'pi_c': 'cpd00009_c0',
    'tagd_c': 'cpd00437_c0',
    'tagd_e' : 'cpd00437_e0',
    'ppat_c' : 'cpd04099_c0',
    'ppat_e' : 'cpd04099_e0',
    

}

conversions_bio_sp = {'acoa_c': 'cpd00022_c0',
 'alaL_c': 'cpd00035_c0',
 'argL_c': 'cpd00051_c0',
 'aspL_c': 'cpd00041_c0', 
 'atp_c': 'cpd00002_c0',
 'citrL_c': 'cpd00274_c0',
 'coa_c': 'cpd00010_c0',
 'fad_c': 'cpd00015_c0',
 'glcAD_c': 'cpd00190_c0',
 'glnL_c': 'cpd00053_c0',
 'gluL_c': 'cpd00023_c0',
 'gly_c': 'cpd00033_c0',
 'hisL_c': 'cpd00119_c0',
 'ileL_c': 'cpd00322_c0',
 'leuL_c': 'cpd00107_c0',
 'lysL_c': 'cpd00039_c0',
 'metL_c': 'cpd00060_c0',
 'nad_c': 'cpd00003_c0',
 'nadh_c': 'cpd00004_c0',
 'nadp_c': 'cpd00006_c0',
 'nadph_c': 'cpd00005_c0',
 'pheL_c': 'cpd00066_c0',
 'pi_c': 'cpd00009_c0',
 'proL_c': 'cpd00129_c0',
 'serL_c': 'cpd00054_c0',
 'thrL_c': 'cpd00161_c0',
 'thym_c': 'cpd00151_c0',
 'valL_c': 'cpd00156_c0',
 'SPs_c': "SPs_c"
}

conversions_bio_lip={
'Lip_c': "Lip_c",
 'dgludmygl_c': "cpd15729_c0",
 'dgludpalmgl_c': "cpd15728_c0",
 'mglucsyldmygl_c': "cpd15738_c0",
 'mglucsyldpalmgl_c': "cpd15737_c0",
 'mycdlpn_c': "cpd15792_c0",
 'myrphgl_c': "cpd15783_c0",
 'palmcdlpn_c': "cpd15791_c0",
 'palmphgl_c': "cpd15782_c0",
 'phosglcdihexdec_c': "cpd15538_c0",
 'phosglcdioctdec_c': "cpd15540_c0",
 'phosglcditetdec_c': "cpd15536_c0",
 'strcdlpn_c': "cpd15793_c0",
 'strphglc_c': "cpd15784_c0"
 }

conversions_bio_prot={
 'Prot_c': 'Prot_c',
  'adp_c': 'cpd00008_c0',
  'alaL_c': 'cpd00035_c0',
  'argL_c': 'cpd00051_c0',
  'asnL_c': 'cpd00132_c0',
  'aspL_c': 'cpd00041_c0',
  'atp_c': 'cpd00002_c0',
  'cysL_c': 'cpd00084_c0',
  'glnL_c': 'cpd00053_c0',
  'gluL_c': 'cpd00023_c0',
  'gly_c': 'cpd00033_c0',
  'h2o_c': 'cpd00001_c0',
  'hisL_c': 'cpd00119_c0',
  'ileL_c': 'cpd00322_c0',
  'leuL_c': 'cpd00107_c0',
  'lysL_c': 'cpd00039_c0',
  'metL_c': 'cpd00060_c0',
  'pheL_c': 'cpd00066_c0',
  'pi_c': 'cpd00009_c0',
  'proL_c': 'cpd00129_c0',
  'serL_c': 'cpd00054_c0',
  'thrL_c': 'cpd00161_c0',
  'trpL_c': 'cpd00065_c0',
  'tyrL_c': 'cpd00069_c0',
  'valL_c': 'cpd00156_c0'
}

conversions_bio_clp={'CLP_c': "CLP_c",
  'adp_c': "cpd00008_c0",
  'alaD_c': "cpd00117_c0",
  'atp_c': "cpd00002_c0",
  'h2o_c': "cpd00001_c0",
  'pi_c': "cpd00009_c0",
  'udp_c': "cpd00014_c0",
  'udpaga_c': "cpd00037_c0",
  'udpamagdapaa_c': "cpd02932_c0"
}

conversions_bio_dna =   {
  'DNA_c': 'DNA_c',
  'adp_c': 'cpd00008_c0',
  'atp_c': 'cpd00002_c0',
  'datp_c': 'cpd00115_c0',
  'dctp_c': 'cpd00356_c0',
  'dgtp_c': 'cpd00241_c0',
  'dttp_c': 'cpd00357_c0',
  'h2o_c': 'cpd00001_c0',
  'pi_c': 'cpd00009_c0',
  'ppi_c': 'cpd00012_c0'
}

conversions_bio_rna =   {'RNA_c': 'RNA_c',
  'adp_c': 'cpd00008_c0',
  'atp_c': 'cpd00002_c0',
  'ctp_c': 'cpd00052_c0',
  'gtp_c': 'cpd00038_c0',
  'h2o_c': 'cpd00001_c0',
  'pi_c': 'cpd00009_c0',
  'ppi_c': 'cpd00012_c0',
  'utp_c': 'cpd00062_c0'
  }

conversions_bio_cw =   {'CLP_c': 'CLP_c',
  'CW_c': "CW_c",
  'teichoic_c': "cpd11442_c0"
}

conversions_bio1 = {'adp_c': 'cpd00008_c0',
  'atp_c': 'cpd00002_c0',
  'CW_c': 'CW_c',
  'DNA_c': 'DNA_c',
  'Lip_c': 'Lip_c',
  'Prot_c': 'Prot_c',
  'RNA_c': 'RNA_c',
  'SPs_c': 'SPs_c',
  'biomass_c': 'biomass_c',
  'h2o_c': 'cpd00001_c0',
  'pi_c': 'cpd00009_c0'
  }

