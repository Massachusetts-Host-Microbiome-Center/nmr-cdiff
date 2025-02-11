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

def parse_stoichiometry(stoich_string):
    stoich_dict = {}
    print(stoich_string)
    for compound in stoich_string.strip("'").split(";"):
        parts = compound.split(":")
        if len(parts) >= 2:
            stoich_dict[parts[1]] = float(parts[0])

    return(stoich_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='bifermentans_model_adding_metabolomics.py')
    parser.add_argument("--params_file", type = str)
    args = parser.parse_args()

    # 0. Read in config file
    config = ConfigParser()
    config.read(args.params_file)
    project_dir = config.get('paths', 'project_dir')


    # 2. Read in model
    model = cobra.io.json.load_json_model(project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff.json")

    # read in modelSEED reaction database
    modelseed_reactions = pd.read_csv("../../data/model_seed_database/ModelSEEDDatabase/Biochemistry/reactions.tsv", sep="\t", index_col = 0)

    # Read in reactions to add
    rxns = pd.read_csv(project_dir + config.get('paths', 'fname_reactions_to_add_metabolomics'))["RXN"]
    rxns = rxns.dropna().str.strip()
    rxns = rxns.values
    rxns = [i for i in rxns if i in modelseed_reactions.index.values]

    # 3. Get the stoichiometry of all reactions to add
    stoichs = [parse_stoichiometry(i) for i in modelseed_reactions.loc[rxns, 'stoichiometry']]

    # 4. Add reactions to the model
    already_present = []
    for rxn, stoich in zip(rxns, stoichs):
        if rxn not in model.reactions:
            new_reaction = Reaction(rxn)
            new_reaction.name = modelseed_reactions.loc[rxn, 'name']
            if pd.isna(new_reaction.name):
                new_reaction.name = "placeholder"
            new_reaction.lower_bound = 0.0  # Assuming default lower bound
            new_reaction.upper_bound = 1000.0  # Assuming default upper bound
            new_reaction.add_metabolites({Metabolite(met_id + "_c0"): coeff for met_id, coeff in stoich.items()})
            model.add_reactions([new_reaction])
        else:
            already_present.append(rxn)
    print("# of reactions in model: ", len(model.reactions))

    print(already_present)

    # 10. Save the model
    cobra.io.json.save_json_model(model.model, project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics.json")
    cobra.io.write_sbml_model(model.model,  project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics.xml")

