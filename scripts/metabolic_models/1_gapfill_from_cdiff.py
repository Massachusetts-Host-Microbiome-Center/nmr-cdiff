import pandas as pd
import numpy as np
#import modelseedpy 
import json
import re
from bioservices import KEGG
import re
from scipy.stats import ranksums
import statsmodels.stats.multitest as multi
import math
from copy import deepcopy
import working_w_seed_models as mm
from importlib import reload
import argparse
import configparser

import cobra
from cobra.medium import minimal_medium
from cobra.io import load_json_model, save_json_model
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import pickle



def find_producing_reactions(model, reaction_id):
    # Get the specified reaction
    reaction = model.reactions.get_by_id(reaction_id)

    # Get the metabolites consumed by the reaction (reactants)
    reactants = reaction.reactants

    # Find reactions that produce these reactants
    producing_reactions = set()
    for reactant in reactants:
        for producing_reaction in reactant.reactions:
            if reactant in producing_reaction.products:
                if producing_reaction.id != reaction_id:  # Exclude the original reaction
                    producing_reactions.add(producing_reaction)
    return producing_reactions


def add_cofactor_synthesis(model, cdiff_model, update_reachable = True):
    # fmnh2 
    reaction_ids = ["ID_692", "ID_692_rev", "ID_263", "ID_170", "ID_170_rev", "ID_460", "ID_460_rev"]
    for reaction_id in reaction_ids:
        model.addCdiffReaction(cdiff_model, reaction_id, enable_exchange=False, check_gene_rule=False, update_reachable = update_reachable)

    # coa
    coa_synth = ['ID_570', 'ID_464', 'ID_139', 'ID_644', 'ID_134'] # from google and cdiff model
    for reaction_id in coa_synth:
        model.addCdiffReaction(cdiff_model, reaction_id, enable_exchange= False, check_gene_rule=False, update_reachable = update_reachable)
    model.addCdiffReaction(cdiff_model, "Trans_cysL", enable_exchange = True, check_gene_rule=False, update_reachable = update_reachable)

    # adding de nove nad+ synthesis
    reaction_ids = ['Trans_nac', 'ID_262', 'ID_27', 'ID_218', 'ID_313']
    for reaction_id in reaction_ids:
        model.addCdiffReaction(cdiff_model, reaction_id, enable_exchange=False, check_gene_rule=False, update_reachable = update_reachable)

    # tetrahydrofolate synthesis
    reaction_ids = ["Trans_chor", "ID_495", "ID_438", "ID_344", "ID_75", "ID_591", "ID_341", "ID_457"]
    # 4-amino-4-deoxychorismate synthase
    # 4-amino-4-deoxychorismate pyruvate-lyase
    # dihydropteroate synthase
    # dihydrofolate synthase
    # "tetrahydrofolate:NADP+ oxidoreductase"
    for reaction_id in reaction_ids:
        model.addCdiffReaction(cdiff_model, reaction_id, enable_exchange=False, check_gene_rule=False, update_reachable = update_reachable)


    # add some limited ferrodoxin import, because c. diff cannot synthesize ferrodoxin to my knowledge
    # TODO: It may be that C. diff and PBI actually can synthesize ferrodoxin, but it's not in the cdiff model
    model.addImportReaction(cdiff_model.conversions['feroxoxi_c'], reaction_id = "Trans_ferr", reaction_name = "Manually transport ferrodoxin",
                        ex_metabolite_id = cdiff_model.conversions['feroxoxi_c'].replace("_c0", "_e0"),
                        ex_metabolite_name = "Ferrodoxin e",
                        in_metabolite_name = "Ferrodoxin c",
                        enable_exchange=True,
                        exchange_lim=2)

    return(model)

def add_biolipid_reactions(model, cdiff_model):
    model.addCdiffReaction(cdiff_model, "ID_05457", enable_exchange=False, check_gene_rule=False, update_reachable=False) #myrstcoa_c
    model.addCdiffReaction(cdiff_model, "ID_05458", enable_exchange=False, check_gene_rule=False, update_reachable=False) #palmcoa_c


    # these were added to cdiff without gene rules by another lab
    model.addCdiffReaction(cdiff_model, "ID_185_1", enable_exchange=False, check_gene_rule=False, update_reachable=False) #phosglcdihexdec_c
    model.addCdiffReaction(cdiff_model, "ID_185_2", enable_exchange=False, check_gene_rule=False, update_reachable=False) #phosglcditetdec_c
    model.addCdiffReaction(cdiff_model, "ID_185_3", enable_exchange=False, check_gene_rule=False, update_reachable=False) #phosglcdioctdec_c

    model.addCdiffReaction(cdiff_model, "ID_200_1", enable_exchange=False, check_gene_rule=False, update_reachable=False) #mpalmphgl_c
    model.addCdiffReaction(cdiff_model, "ID_200_1_rev", enable_exchange=False, check_gene_rule=False, update_reachable=False) #mpalmphgl_c
    model.addCdiffReaction(cdiff_model, "ID_200_2", enable_exchange=False, check_gene_rule=False, update_reachable=False) #myrphgl_c
    model.addCdiffReaction(cdiff_model, "ID_200_2_rev", enable_exchange=False, check_gene_rule=False, update_reachable=False) #myrphgl_c
    model.addCdiffReaction(cdiff_model, "ID_200_3", enable_exchange=False, check_gene_rule=False, update_reachable=False) #strphglc_c
    model.addCdiffReaction(cdiff_model, "ID_200_3_rev", enable_exchange=False, check_gene_rule=False, update_reachable=True) #strphglc_c
    model.addCdiffReaction(cdiff_model, "ID_9624", enable_exchange=False, check_gene_rule=False, update_reachable=True) #srcoa
    model.addCdiffReaction(cdiff_model, "ID_9624_rev", enable_exchange=False, check_gene_rule=False, update_reachable=True) #srcoa
    return(model)


def check_biomass_reactions(model, cdiff_model):
    model.enableTransportReactions()
    model.model.objective = "Bio_DNA"
    model.addSink(cdiff_model.conversions["DNA_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce DNA!")
    else:
        print("Model IS NOT able to produce DNA")


    model.enableTransportReactions()
    model.model.objective = "Bio_RNA"
    model.addSink(cdiff_model.conversions["RNA_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce RNA!")
    else:
        print("Model IS NOT able to produce RNA")


    model.model.objective = "Bio_CW"
    model.addSink(cdiff_model.conversions["CW_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce CW!")
    else:
        print("Model IS NOT able to produce CW")

    model.model.objective = "Bio_prot"
    model.addSink(cdiff_model.conversions["Prot_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce Protein!")
    else:
        print("Model IS NOT able to produce Protein")

    model.model.objective = "Bio_lip"
    model.addSink(cdiff_model.conversions["Lip_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce lipid!")
    else:
        print("Model IS NOT able to produce lipid")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='bifermentans_model.py')
    parser.add_argument("--params_file", type = str)
    parser.add_argument("--build_model", action = 'store_true')
    args = parser.parse_args()

    # 0. Read in config file
    config = configparser.ConfigParser()
    config.read(args.params_file)
    project_dir = config.get('paths', 'project_dir')
    
    if args.build_model:
        # 1. Load both PBI and C. diff model
        model_cobra_cdiff = cobra.io.load_json_model(project_dir + config.get('paths', 'model_cdiff'))
        cdiff_model = mm.Model(model_cobra_cdiff, config) # automatically loads conversions into cdiff_model.conversions

        model_cobra = cobra.io.load_json_model(project_dir + config.get('paths', 'model_query'))
        model = mm.Model(model_cobra, config)
        print("Number of reactions in raw: ", len(model.model.reactions))
        print("Number of metabolites in raw: ", len(model.model.metabolites))

        # 2. Read in the gene rule conversion from PBI to C. diff based on the cell host microbe paper, and add those gene rules to the model object (used if check_gene_rule = True)
        pbi_cdiff_generule = pd.read_excel(project_dir + config.get('paths', 'gene_rule_pbi_cdiff'), sheet_name = "GC.4 PBI Gene Map")
        model.gene_rule_df = pbi_cdiff_generule

        # 3. Add relevant biomass reactions from cdiff model
        model.addCdiffReaction(cdiff_model, 'Bio_CW', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_prot', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_CLP', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_RNA', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_DNA', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_SPs', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'Bio_lip', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.addCdiffReaction(cdiff_model, 'bio1', enable_exchange = False, check_gene_rule=False, update_reachable = False)
        model.model.objective = "bio1"

        print("Number of reactions in adding biomass endpoints: ", len(model.model.reactions))
        print("Number of metabolites in adding biomass endpoints: ", len(model.model.metabolites))


        # 4. Add reactions from cdiff model to PBI model based on gene rules
        gene_rules_pbi = model.gene_rule_df['CD630_homolog'].values
        for reaction in cdiff_model.model.reactions:
            # single gene rule
            gene_rule = reaction.gene_reaction_rule
            if gene_rule in gene_rules_pbi:
                model.addCdiffReaction(cdiff_model, reaction.id, check_gene_rule = True, update_reachable = False)
                
            # stacked gene _rules
            gene_rules = []
            if 'or' in gene_rule:
                gene_rules = gene_rule.split(" or ")
            if 'and' in gene_rule:
                gene_rules = gene_rule.split(" and ")
            for gene_rule in gene_rules:
                if gene_rule in gene_rules_pbi:
                    model.addCdiffReaction(cdiff_model, reaction.id, check_gene_rule = True, update_reachable = False)

        print("Number of reactions in adding gene rule reactions: ", len(model.model.reactions))
        print("Number of metabolites in adding gene rule reactions: ", len(model.model.metabolites))


        # 5. Calculate reachability

        transport_reactions = model.getTransportReactions()
        model.enableTransportReactions(transport_reactions)
        previously_reachable_metabolites = -1
        while len(model.reachable_metabolite_ids) != previously_reachable_metabolites:
            previously_reachable_metabolites = len(model.reachable_metabolite_ids)
            model.reachable_reaction_ids, not_reachable_reaction_ids, model.reachable_metabolite_ids, not_reachable_metabolite_ids = model.findUnreachables()

        print("Number of reactions reachable: ", len(model.reachable_reaction_ids))
        print("% reachable reactions: ", len(model.reachable_reaction_ids) / len(model.model.reactions))


        # 6. add reactions to model for synthesizing cofactors. The cofactors are missing from the model, and we know that from C. diff model that they are likely to be synthesized by these organisms
        model = add_cofactor_synthesis(model, cdiff_model)
        model.reachable_reaction_ids, not_reachable_reaction_ids, model.reachable_metabolite_ids, not_reachable_metabolite_ids = model.findUnreachables()
        print("Number of reactions in adding cofactor synthesis: ", len(model.model.reactions))
        print("Number of reactions reachable: ", len(model.reachable_reaction_ids))
        print("% reachable reactions: ", len(model.reachable_reaction_ids) / len(model.model.reactions))

        # 7. manually add reactions necessary to reach Bio_CW. Sometimes the tree branching strategy hits an infinite loop, so we manually add reactions necessary to reach Bio_CW.
        reaction_ids = ["ID_626", "ID_462", "ID_59", "ID_607", "ID_324", 'Trans_cdpg']
        for reaction_id in reaction_ids:
            model.addCdiffReaction(cdiff_model, reaction_id, enable_exchange=False, check_gene_rule=False, update_reachable=True)
        model.reachable_reaction_ids, not_reachable_reaction_ids, model.reachable_metabolite_ids, not_reachable_metabolite_ids = model.findUnreachables()
        print("Number of reactions in adding BioCW manuals: ", len(model.model.reactions))
        print("Number of reactions reachable: ", len(model.reachable_reaction_ids))
        print("% reachable reactions: ", len(model.reachable_reaction_ids) / len(model.model.reactions))

        # 8. Add fatty acid import reactions
        model.addImportReaction(cdiff_model.conversions['myrstacp_c'], reaction_id = "Trans_myrst_acp", reaction_name = "Manually transport myrstic acid acp",
                        ex_metabolite_id = cdiff_model.conversions['myrstacp_c'].replace("_c0", "_e0"),
                        ex_metabolite_name = "Myristic acid acp e",
                        in_metabolite_name = "Myristic acid acp c",
                        enable_exchange=True,
                        exchange_lim=10)

        # Add fatty acid transporters - this may have to change to be fatty acids directly and not CoA
        model.addImportReaction(cdiff_model.conversions['palmacp_c'], reaction_id = "Trans_palm_acp", reaction_name = "Manually transport palmitic acid acp",
                            ex_metabolite_id = cdiff_model.conversions['palmacp_c'].replace("_c0", "_e0"),
                            ex_metabolite_name = "Palmitic acid acp e",
                            in_metabolite_name = "Palmitic acid acp c",
                            enable_exchange=True,
                            exchange_lim=10)

        # Add fatty acid transporters - this may have to change to be fatty acids directly and not CoA
        model.addImportReaction(cdiff_model.conversions['sracp_c'], reaction_id = "Trans_sr_acp", reaction_name = "Manually transport Stearoyl acp",
                            ex_metabolite_id = cdiff_model.conversions['sracp_c'].replace("_c0", "_e0"),
                            ex_metabolite_name = "Stearoyl acp e",
                            in_metabolite_name = "Stearoyl acp c",
                            enable_exchange=True,
                            exchange_lim=10)
        
        print("Number of reactions in adding fatty acid import reactions: ", len(model.model.reactions))
        print("Number of reactions reachable: ", len(model.reachable_reaction_ids))
        print("% reachable reactions: ", len(model.reachable_reaction_ids) / len(model.model.reactions))

        # 9. Add biolipid reactions
        model = add_biolipid_reactions(model, cdiff_model)

        # 10. Add reactions using tree branching strategy
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_CW", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_prot", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_CLP", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_SPs", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_lip", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_RNA", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "Bio_DNA", check_gene_rule=False, cycle_lim = 10)
        reactions_added = model.gapfill_from_cdiff_tree(cdiff_model, "bio1", check_gene_rule=False, cycle_lim = 10)
        

        model.reachable_reaction_ids, not_reachable_reaction_ids, model.reachable_metabolite_ids, not_reachable_metabolite_ids = model.findUnreachables()
        print("Number of reactions in adding BioCW manuals: ", len(model.model.reactions))
        print("Number of reactions reachable: ", len(model.reachable_reaction_ids))
        print("% reachable reactions: ", len(model.reachable_reaction_ids) / len(model.model.reactions))

        
        check_biomass_reactions(model, cdiff_model)
        

        model.model.objective = "bio1"
        model.addSink(cdiff_model.conversions["biomass_c"])
        if model.model.optimize().objective_value > 0:
            print("\n Model is able to produce Biomass!")
        else:
            print("Model IS NOT able to produce Biomass")

        pickle.dump(model, open(project_dir + config.get('paths', 'output_dir') + "bifermentans_gapfilled_from_cdiff.pkl", "wb"))
        cobra.io.save_json_model(model.model, project_dir + config.get('paths', 'output_dir') + "bifermentans_gapfilled_from_cdiff.json") 


    #########################################################
    ############## PRINT OUT MODEL CHANGES  #################
    #########################################################

    def addGeneRuleAnnotation(df, cdiff_model, config):
        filepath = project_dir + config.get('paths', 'gene_rule_pbi_cdiff')
        gene_rule_df = pd.read_excel(filepath, sheet_name='GC.4 PBI Gene Map')
        gene_rule_df.index = gene_rule_df['CD630_homolog']
        gene_rule_df = gene_rule_df.dropna(subset=["CD630_homolog"])
        gene_rule_df = gene_rule_df.drop_duplicates(subset=["CD630_homolog"])
        gene_rule_df 

        gene_rules  = [cdiff_model.model.reactions.get_by_id(i.id).gene_reaction_rule if i.id in cdiff_model.reaction_id_list else "" for i in reactions_added]
        gene_rules

        gene_rules = [re.split(' or | and ', i) for i in gene_rules]
        if ' and ' in gene_rules:
            gene_rules = [i.split(' and ') for i in gene_rules]


        for gene_rule_list, reaction in zip(gene_rules, reactions_added):
            for g in gene_rule_list:
                if g in gene_rule_df.index.values:
                        df.loc[reaction.id, "CD630_annotation"] = gene_rule_df.loc[g, "CD630_homolog_annotation"]
                        df.loc[reaction.id, "PBI_annotation"] = gene_rule_df.loc[g, "PATRIC_annotation"]
        return(df)
    
    def getReactionString(model, reaction_id):
        reaction = model.reactions.get_by_id(reaction_id)
        reactant_names = [i.name.replace(" [c0]", "") for i in reaction.reactants]
        product_names = [i.name.replace(" [c0]", "") for i in reaction.products]
        reaction_str = ' + '.join(reactant_names) + "=>" + " + ".join(product_names)
        return(reaction_str)

    model = pickle.load(open(project_dir + config.get('paths', 'output_dir') + "bifermentans_gapfilled_from_cdiff.pkl", "rb"))
    model_gapfill = model.model
    model_orig = cobra.io.load_json_model(project_dir + config.get('paths', 'model_query'))
    model_cobra_cdiff = cobra.io.load_json_model(project_dir + config.get('paths', 'model_cdiff'))
    cdiff_model = mm.Model(model_cobra_cdiff, config) # automatically loads conversions into cdiff_model.conversions

    check_biomass_reactions(model, cdiff_model)
    model.model.objective = "bio1"
    model.addSink(model.conversions["biomass_c"])
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce Biomass!")
    else:
        print("Model IS NOT able to produce Biomass")

    reactions_added = [i for i in model_gapfill.reactions if i not in model_orig.reactions]
    reactions_added = [i for i in reactions_added if "SK_" not in i.id and "test" not in i.id]
    reactions_added = [i for i in reactions_added if "EX_" not in i.id]
    ids = [i.id for i in reactions_added]
    names = [i.name for i in reactions_added]
    reaction_strings = [getReactionString(model_gapfill, i.id) for i in reactions_added]
    cdiff_model.reaction_id_list = [i.id for i in cdiff_model.model.reactions]
    gene_rules  = [cdiff_model.model.reactions.get_by_id(i.id).gene_reaction_rule if i.id in cdiff_model.reaction_id_list else "" for i in reactions_added]

    df = pd.DataFrame({'id': ids, 'name': names, 'reaction_string': reaction_strings, 'gene_rule': gene_rules})
    df.index = df['id']
    #df.to_csv(project_dir + config.get('paths', 'output_dir') + "bifermentans_gapfilled_from_cdiff_changes.csv")

    df = addGeneRuleAnnotation(df, cdiff_model, config)
    df.to_csv(f"{project_dir}/{config.get('paths', 'output_dir')}/reactions_added_gapfill.csv")


