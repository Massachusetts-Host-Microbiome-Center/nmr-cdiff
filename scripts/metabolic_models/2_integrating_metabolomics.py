import pandas as pd
import numpy as np
import modelseedpy 
import json
import re
import re
import cobra
from cobra.io import load_json_model, save_json_model
from cobra import Model, Reaction, Metabolite
from bioservices import KEGG
import re
from scipy.stats import ranksums
import statsmodels.stats.multitest as multi
import math
from copy import deepcopy
import scipy
import seaborn as sns

import working_w_seed_models as mm
from importlib import reload
from cobra.flux_analysis import flux_variability_analysis
import configparser
import pickle
import argparse
from matplotlib import pyplot as plt



# manual metabolomic to model seed conversions
d2 = {'2-hydroxybutyrate/2-hydroxyisobutyrate':'cpd03561',
      '4-hydroxyphenylacetate sulfate\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0':'cpd00489',
      '5-methylthioribose**':'cpd01981',
      'bilirubin (Z,Z)':'cpd00376',
      'cis-urocanate':'cpd00581',
      'riboflavin (Vitamin B2)': 'cpd00220',
      'tauroursodeoxycholate': 'cpd17168',
      'pipecolate': 'cpd00323',
      'sphingadienine': 'cpd25264',
      'palmitoylcarnitine': 'cpd01915',
      'oleoylcarnitine (C18:1)': 'cpd33454',
      'nicotinate ribonucleoside': 'cpd03471',
      'N-formylmethionine': 'cpd02012',
      "N('1)-acetylspermidine": 'cpd00470',
      'methylmalonate (MMA)':'cpd01468',
      'gamma-aminobutyrate (GABA)':'cpd11408',
       'hexadecasphingosine (d16:1)*':'',
       'laurate (12:0)':'cpd01741',
      'N-acetylglutamine': 'cpd01760',
      'N-acetyl-beta-glucosaminylamine': 'cpd00912',
     'N-acetylasparagine':'cpd25519', 
     'N-acetylaspartate (NAA)':'cpd00767',
     'N-acetylcitrulline': 'cpd11209',
     'N-acetylglucosaminylasparagine': 'cpd02763',
     'N-acetylglutamate': 'cpd00477',
     'N-acetylhistidine': 'cpd24478',
     'N-acetyltryptophan': 'cpd27567',
     'threonine': 'cpd00161',
     'proline': 'cpd00129',
     'N1-methylguanosine':'cpd25947',
     'N1-methyladenosine': 'cpd01637',
     'histidine betaine (hercynine)*': 'cpd03305',
     'glycerophosphoserine*': 'cpd15468',
     'cysteine': 'cpd00084'

     }


# Read in metabolomics data
def getData(fname_metabolomics):
    df = pd.read_excel(fname_metabolomics, sheet_name="GF.5 MetabolomicsDataset", skiprows = np.arange(7), index_col = "BIOCHEMICAL")

    feature_info = ["PATHWAY_SORTORDER", "SUPER_PATHWAY", "SUB_PATHWAY", "COMP_ID", "PLATFORM", "CHEMICAL_ID", "RI", "MASS", "CAS", "PUBCHEM", "CHEMSPIDER", "KEGG","                               Group   HMDB_ID"]
    feature_annotation = df.loc[:, feature_info]
    df = df.drop(feature_info, axis = 1)

    # Drop rows (metabolites) with too many NAs
    print("Starting with ", df.shape[0], " features")
    keep = df.isna().sum(axis = 1) / df.shape[1] < 0.5
    print("Keeping ", sum(keep), " after removing those with > 25% NA")
    df = df.loc[keep, :]

    # Drop columns (samples) with too many NAs
    print("Starting with ", df.shape[1], " samples")
    keep = df.isna().sum(axis = 0) / df.shape[0] < 0.25
    print("Keeping ", sum(keep), " after removing those with > 25% NA")
    df = df.loc[:, keep]
    
    
    df_norm = df.div(df.median(axis=1), axis=0)
    df_norm = df_norm.transform('log')
    df_norm = df_norm.loc[["X -" not in i for i in df_norm.index.values], :]
    
    df = df.transform('log')
    df = df.loc[["X -" not in i for i in df.index.values], :]
    
    print(df.shape)
    return(df.T, df_norm.T)

# read in metadata

def getMetadata(df):
    metadata = pd.DataFrame(index = df.index.values)
    metadata['condition'] = [str(i).split(" ")[0] for i in df.index.values]
    metadata['condition'] = [str(i).split(".")[0] for i in metadata['condition'].values]
    metadata['condition'][metadata['condition'] == "Germ"] = "GF"
    metadata.condition.unique()
    return(metadata)

def findDifferentialMetabolites(df_norm, condition1, condition2, alternative):

	gf = df_norm.loc[df_norm.index.values == condition1, ]
	pbi = df_norm.loc[df_norm.index.values == condition2, ]
	# calc. metabolites with a significant increase in pbi

	pvals = []
	for met in pbi.columns.values:
	    tmp = scipy.stats.ranksums(gf.loc[:, met].values, pbi.loc[:, met].values, alternative = alternative)
	    pvals.append(tmp[1])
	pvals = np.array(pvals)
	pvals[[math.isnan(i) for i in pvals]] = 1

	p_adj = multi.multipletests(pvals, method = "fdr_by")[1]
	sns.histplot(p_adj, bins = 100)
	mets = pbi.columns.values[p_adj <= .05]

	return(mets)

def metaboliteIsReachable(model, met):
    # finds all reactions that produce the metabolite
    # removes any direct transport reactions that produce the metabolite
    # checks if metabolite is still reachable
    
    producing_reactions = model.find_producing_reactions(met)
    direct_transport = [i for i in producing_reactions if i.id in model.getTransportReactions()]
    if len(direct_transport) > 0:
        reaction_save = direct_transport
        model.model.remove_reactions(direct_transport)
    producing_reactions = model.find_producing_reactions(met)
    transport_reactions = model.getTransportReactions()
    model.enableTransportReactions()
    reachable_reaction_ids, not_reachable_reactions, reachable_metabolite_ids, not_reachable_metabolites = model.findUnreachables(transport_reactions)
    return(met in reachable_metabolite_ids)


def calculateReachableReactions(model):
	# calc reachable
	num_reachable_metabolites = []
	num_reachable_reactions = []
	reachable_metabolites = {}
	previously_reachable_metabolites = -1
	i = 0
	while len(model.reachable_metabolite_ids) != previously_reachable_metabolites:
	    previously_reachable_metabolites = len(model.reachable_metabolite_ids)
	    transport_reactions = model.getTransportReactions()
	    model.enableTransportReactions()
	    reachable_reaction_ids, not_reachable_reactions, reachable_metabolite_ids, not_reachable_metabolites = model.findUnreachables(transport_reactions)
	    model.reachable_metabolite_ids = reachable_metabolite_ids
	    model.reachable_reaction_ids = reachable_reaction_ids
	    model.reachable_metabolite_ids
	    num_reachable_metabolites.append(len(model.reachable_metabolite_ids))
	    num_reachable_reactions.append(len(model.reachable_reaction_ids))
	    reachable_metabolites[i] = model.reachable_metabolite_ids
	    print(i, len(model.reachable_metabolite_ids))
	    i = i + 1
	return(model)

def convertMetabolomicsToModelSeed(mets, fname_compound_lookup):
	compound_df = pd.read_csv(fname_compound_lookup, sep = "\t")
	keep = [i == i for i in compound_df['aliases']]
	compound_df = compound_df.loc[keep, :]

	d = {}
	#d2 = {}
	for met in mets:
	    for i in range(compound_df.shape[0]):
	        aliases = compound_df.aliases.values[i].split(";")
	        aliases = [i.replace("Name: ", "").strip().lower() for i in aliases]
	        if met in aliases:
	            d[met] = compound_df.id.values[i]
	            #d2[met] = aliases
	return(d)

def plot_metabolite_differential(df_norm, mets, condition1, condition2):
	df_norm = pd.concat([df_norm.loc[df_norm.index.values == condition1, ], df_norm.loc[df_norm.index.values == condition2, ]]).loc[:, mets]
	df_norm['condition'] = df_norm.index.values
	df_melt = df_norm.melt('condition')
	df_melt['BIOCHEMICAL'] = [i.replace("BIOCHEMICAL=", "") for i in df_melt.BIOCHEMICAL.values]

	#sns.boxplot(df, facet_kws=dict(margin_titles=True))
	g = sns.FacetGrid(df_melt, col="BIOCHEMICAL", sharey= False, col_wrap = 5)
	g.map(sns.boxplot, "condition", "value")
	g.set_titles("{col_name}");

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='bifermentans_model_adding_metabolomics.py')
	parser.add_argument("--params_file", type = str)
	args = parser.parse_args()

    # 0. Read in config file
	config = configparser.ConfigParser()
	config.read(args.params_file)
	project_dir = config.get('paths', 'project_dir')

	# 1. Read in data
	params_file = args.params_file
	config = configparser.ConfigParser()
	print(params_file)
	config.read(params_file)
	model_cobra = cobra.io.load_json_model(project_dir + config.get('paths', 'output_dir') + "bifermentans_gapfilled_from_cdiff.json")
	model = mm.Model(model_cobra, config)
	model.model.optimize().objective_value

	df, df_norm = getData(project_dir + config.get('paths', 'fname_metabolomics'))
	metadata = getMetadata(df)
	df.index = metadata.loc[df.index.values, "condition"]
	df_norm.index = metadata.loc[df_norm.index.values, "condition"]

	condition1 = "GF"
	condition2 = "PBI"		
	mets_used = findDifferentialMetabolites(df_norm, condition1, condition2, "greater")
	mets_produced = findDifferentialMetabolites(df_norm, condition1, condition2, "less")


	plot = True
	if plot:
		p = plot_metabolite_differential(df_norm, mets_produced, condition1, condition2)
		plt.savefig(project_dir + config.get('paths', 'output_dir') + "/differential_metabolites_produced.png")
		p = plot_metabolite_differential(df_norm, mets_used, condition1, condition2)
		plt.savefig(project_dir + config.get('paths', 'output_dir') + "/differential_metabolites_utilized.png")


	d = convertMetabolomicsToModelSeed(mets_produced, project_dir + config.get('paths', 'modelseed_compound_db'))
	met_id_dict = {**d, **d2}

	print("Number metabolites increased in vitro but missing modelseed ID: ", len([i for i in mets_produced if i not in met_id_dict]))
	print("Number of metabolites increased in vitro and have modelseed ID: ", len([i for i in mets_produced if i in met_id_dict]))
	print("Missing ids: ")
	[i for i in mets_produced if i not in met_id_dict]

	# 1. For compounds already in the model, add reactions from c. diff model to make them reachable
	## A. Identify compounds that are already in the model=
	need_to_reach = [met_id_dict[i] + "_c0" for i in mets_produced if i in met_id_dict]
	need_to_reach = np.array(need_to_reach)[[i in model.metabolite_id_list for i in need_to_reach]]
	print("Metabolites that are produced in metabolomics and already in the model: ", [model.model.metabolites.get_by_id(i).name for i in need_to_reach])

	## B. Identify which of those compounds are already reachable
	model = calculateReachableReactions(model)

	accessible = [(i, model.model.metabolites.get_by_id(i).name) for i in need_to_reach if metaboliteIsReachable(model, i)]
	accessible = pd.DataFrame(accessible, columns=["met_id", "name"])
	print("Metabolites already reachable: ", accessible)

	missing = [(i, model.model.metabolites.get_by_id(i).name) for i in need_to_reach if not metaboliteIsReachable(model, i)]
	missing = pd.DataFrame(missing, columns=["met_id", "name"])
	print("Metabolites that need more reactions: ", missing)

	## C. For the metabolites that need more reactions, port from C. diff model where possible
	model_cobra = cobra.io.load_json_model(project_dir + config.get('paths', 'model_cdiff'))
	cdiff_model = mm.Model(model_cobra, config)

	conversions2 = {value: key for key, value in cdiff_model.conversions.items()}
	for met_id in missing.met_id.values:
		if met_id in conversions2:
			producing_reactions = cdiff_model.find_producing_reactions(conversions2[met_id])
			producing_reactions = [i for i in producing_reactions if i.id not in cdiff_model.getTransportReactions()]
			print([i.name for i in producing_reactions])
			model.addCdiffReaction(cdiff_model, producing_reactions[0].id, enable_exchange = False, check_gene_rule=False, update_reachable = True)
			print(producing_reactions[0].id)
			#model.gapfill_from_cdiff_tree(cdiff_model, producing_reactions[0].id, check_gene_rule = False)

	## D. ID metabolites still not reachable after porting
	missing = [(i, model.model.metabolites.get_by_id(i).name) for i in need_to_reach if not metaboliteIsReachable(model, i)]
	missing = pd.DataFrame(missing, columns=["met_id", "name"])
	missing

	## E. Write out metabolites that are not reachable, either because we didn't get there with C. diff reactions, or because they aren't in the model (but they do have an id)
	missing2 = [(met_id_dict[i], i) for i in mets_produced if i in met_id_dict and met_id_dict[i]+"_c0" not in model.metabolite_id_list and "cpd" in met_id_dict[i]+"_c0"]
	missing2 = pd.DataFrame(missing2, columns=["met_id", "name"])
	pd.concat([missing, missing2]).to_csv(project_dir + config.get('paths', 'outdir') + "/compounds_to_be_added_metabolomics.csv")
	print("Metabolites to be added: ", pd.concat([missing, missing2]))


	# F. Write out metabolites that are utilized in metabolomics but not in model:
	d = convertMetabolomicsToModelSeed(mets_used, project_dir + config.get('paths', 'modelseed_compound_db'))
	df = pd.DataFrame.from_dict(d, orient='index').reset_index()
	df.columns = ['metabolite_name', 'compound_id']
	df['compound_id'] = df['compound_id'] + "_c0"
	df['in_model'] = df['compound_id'].isin(model.metabolite_id_list)

	print("Number metabolites increased in vitro but missing modelseed ID: ", len([i for i in mets_produced if i not in met_id_dict]))
	print("Number of metabolites increased in vitro and have modelseed ID: ", len([i for i in mets_produced if i in met_id_dict]))
	print("Missing ids: ")
	[i for i in mets_produced if i not in met_id_dict]