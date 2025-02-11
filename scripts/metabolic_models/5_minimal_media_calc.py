import pandas as pd
import numpy as np
import cobra
from configparser import ConfigParser
import working_w_seed_models as wm
from importlib import reload
import argparse
import re
reload(wm)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='bifermentans_model_adding_metabolomics.py')
    parser.add_argument("--params_file", type = str)
    args = parser.parse_args()

    # 0. Read in config file
    config = ConfigParser()
    config.read(args.params_file)
    project_dir = config.get('paths', 'project_dir')

    config = ConfigParser()
    config.read("build_pbi_model.ini")
    tmp = cobra.io.load_json_model(project_dir + config.get('paths', 'output_dir') + "/bifermentans_gapfilled_from_cdiff_metabolomics_qc.json")
    model = wm.Model(tmp, config)
    model.model.objective = "bio1"

    model.model.objective = "bio1"
    if model.model.optimize().objective_value > 0:
        print("\n Model is able to produce Biomass!")
    else:
        print("Model IS NOT able to produce Biomass")


    # 1. remove double amino acids as options - we will be adding single amino acids

    # Define the pattern to match metabolites with names like "Ala-His [e0]"
    pattern1 = r'^[a-z]{3}-[A-Z]-[a-z]{3}-[A-Z] \[[a-z]\d\] exchange$'
    pattern2 = r'^[A-Z]{1}[a-z]{2}-[A-Z]{1}[a-z]{2} \[[a-z]\d\] exchange$'
    pattern3 = r'^[a-z]{3}-[a-z]{3}-[A-Z] \[[a-z]\d\] exchange$'
    pattern4 = r'^[a-z]{3}-L-[A-Z][a-z]{2}-L \[[a-z]\d\] exchange$'

    # Iterate through the boundary metabolites and remove those that match the pattern
    for reaction in model.model.boundary:
        if re.match(pattern1, reaction.name) or re.match(pattern2, reaction.name) or re.match(pattern3, reaction.name) or re.match(pattern4, reaction.name):
            model.model.remove_reactions([reaction.id])


    # 2. Get ranked list of metabolites based on how often they are included in the minimal media formulation
    min_media_counts = {}
    concentrations = list(np.arange(0.9, 1, 0.005))
    for concentration in concentrations:
        met_ids = model.calculateMinimalMedia(concentration, minimize_components = True) # minimize_components = True means fewest number of reactions (L1). False will minimize total flux (L2)
        mets = [model.model.metabolites.get_by_id(i.replace("EX_", "")) for i in met_ids]
        min_media = [i.name for i in mets]
        
        for metabolite in mets:
            if metabolite in min_media_counts:
                min_media_counts[metabolite] += 1
            else:
                min_media_counts[metabolite] = 1

    # Convert the dictionary into a pandas DataFrame
    min_media_df = pd.DataFrame(list(min_media_counts.items()), columns=['Metabolite ID', 'Count'])

    # Add a column for the metabolite name
    min_media_df['Metabolite Name'] = min_media_df['Metabolite ID'].apply(lambda x: model.model.metabolites.get_by_id(x.id).name)
    min_media_df['Metabolite Formula'] = min_media_df['Metabolite ID'].apply(lambda x: model.model.metabolites.get_by_id(x.id).formula)

    min_media_df.sort_values(by=["Count", "Metabolite Name"], ascending=[False, True], inplace=True)


    # 3. Calculate the flux range (assume glucose is at 30mM)
    # This assumption happens within calculateMinimalMedia
    min_media_df.index = ["EX_" + str(i) for i in min_media_df["Metabolite ID"].values]
    for reaction_id in min_media_df.index:
        model.enableTransportReaction(reaction_id, concentration = 30) # set max concentration of metabolites to 30 mM
    flux_df = model.model.optimize().fluxes
    min_media_df['flux (mM) in 1 hr'] = -flux_df.loc[min_media_df.index]
    min_media_df.to_csv(project_dir + config.get('paths', 'output_dir') + "minimal_media_predictions.csv")

