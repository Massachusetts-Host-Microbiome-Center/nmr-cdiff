import pandas as pd
import numpy as np
import cobra
import os
import json
from cobra import Model, Reaction, Metabolite
import argparse
from configparser import ConfigParser
import working_w_seed_models as wm
from cobra.core.formula import Formula


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='bifermentans_model_adding qc_reactions.py')
    parser.add_argument("--params_file", type = str)
    args = parser.parse_args()

    # 0. Read in config file
    config = ConfigParser()
    config.read(args.params_file)
    project_dir = config.get('paths', 'project_dir')

    # 1.set up models
    tmp = cobra.io.load_json_model(project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics.json")
    model = wm.Model(tmp, config)
    tmp = cobra.io.load_json_model("../../data/cdiff/icdf843_2_unidirectional.json")
    cdiff_model = wm.Model(tmp, config)

    # 2. Add extra reactions discovered during QC
    model.addCdiffReaction(cdiff_model, "Trans_cysL", check_gene_rule=False, update_reachable=False )
    model.addCdiffReaction(cdiff_model, 'Trans_proL', check_gene_rule=False, update_reachable=False)
    model.addCdiffReaction(cdiff_model, 'Trans_proL_PMF', check_gene_rule=False, update_reachable=False)
    model.addCdiffReaction(cdiff_model, "ID_314", check_gene_rule=False, update_reachable=False)
    had_reactions = ['ID_28', 'ID_28_rev', 'ID_685', 'ID_605', 'ID_382', 'ID_31', "ID_20"]
    for reaction in had_reactions :
        model.addCdiffReaction(cdiff_model, reaction, check_gene_rule=False, update_reachable=False)

    # 3. Make cellobiose export 1 directional 
    model.model.reactions.get_by_id("EX_cpd00158_e0").lower_bound = 0 # no import of cellobiose
    model.model.reactions.get_by_id("EX_cpd00324_e0").lower_bound = 0 # no import of mercaptomethane

    # 4. Directionality fixes found during memote qc - doublechecked with modelSEED database
    model.model.reactions.get_by_id("rxn10344_c0").lower_bound = -1000

    # 5. Rename a couple of important metabolites for clarity
    # rename for clarity of final result
    model.model.metabolites.get_by_id("cpd00158_e0").name = "cellobiose [e0]"
    model.model.metabolites.get_by_id("cpd00276_e0").name = "D-Glucosamine [e0]"
    model.model.metabolites.get_by_id("cpd00644_e0").name = "Pantothenic acid [e0]"
    model.model.metabolites.get_by_id("cpd00324_e0").name = "mercaptomethane [e0]"
    

    # 6. Add formulas to metabolites
    compounds = pd.read_csv(project_dir + config.get('paths', 'modelseed_compound_db'), sep="\t", index_col=0)
    compounds.index = [str(i)+"_c0" for i in compounds.index]
    compounds['formula']

    for metabolite in model.model.metabolites:
        if metabolite.id in compounds.index:
            form = compounds.loc[metabolite.id, 'formula']
            if not pd.isna(form):
                metabolite.formula = form
            else:
                metabolite.formula = None
    compounds.index = compounds.index.str.replace("_c0", "_e0")
    for metabolite in model.model.metabolites:
        if metabolite.id in compounds.index:
            form = compounds.loc[metabolite.id, 'formula']
            if not pd.isna(form):
                metabolite.formula = form
            else:
                metabolite.formula = None


    # 7. Remove reactions that use oxygen
    # Find reactions that use oxygen
    reactions_using_cpd00007 = [reaction for reaction in model.model.reactions if 'cpd00007_c0' in [i.id for i in reaction.metabolites]]
    for reaction in reactions_using_cpd00007:
        print(f"Reaction ID: {reaction.id}, Name: {reaction.name}, Metabolites: {reaction.metabolites}")

        model.model.remove_reactions([reaction.id for reaction in reactions_using_cpd00007])


    # 8. Save the model
    cobra.io.json.save_json_model(model.model, project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics_qc.json")
    cobra.io.write_sbml_model(model.model,  project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics_qc.xml")