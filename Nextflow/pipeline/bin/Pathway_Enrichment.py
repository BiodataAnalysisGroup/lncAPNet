#!/usr/bin/env python

import gseapy as gp
import pandas as pd
import numpy as np

import os
import openpyxl
from sklearn.model_selection import train_test_split

import gc

working_dir = "/mnt/c/Users/vvasileiou/Desktop/lncAPNet/"
os.chdir(working_dir)

### PATHWAY ENRICHMENT

# Load the Excel files
icgc = pd.read_excel("Pilot_1_CLL/results/ICGC/ICGC_NetBID2_ms_tab.xlsx")
bcmo = pd.read_excel("Pilot_1_CLL/results/BCMO/BCMO_NetBID2_ms_tab.xlsx")

# Filtering based on conditions
icgc_pos = icgc[(icgc["adj.P.Val.U.Vs.M_DA"] < 0.05) & (icgc["logFC.U.Vs.M_DA"] > 0) & (icgc["Size"] > 30)]
bcmo_pos = bcmo[(bcmo["adj.P.Val.U.Vs.M_DA"] < 0.05) & (bcmo["logFC.U.Vs.M_DA"] > 0) & (bcmo["Size"] > 30)]
icgc_neg = icgc[(icgc["adj.P.Val.U.Vs.M_DA"] < 0.05) & (icgc["logFC.U.Vs.M_DA"] < 0) & (icgc["Size"] > 30)]
bcmo_neg = bcmo[(bcmo["adj.P.Val.U.Vs.M_DA"] < 0.05) & (bcmo["logFC.U.Vs.M_DA"] < 0) & (bcmo["Size"] > 30)]

# Find matching originalID values
pos = icgc_pos[icgc_pos["originalID"].isin(bcmo_pos["originalID"])]
neg = icgc_neg[icgc_neg["originalID"].isin(bcmo_neg["originalID"])]

# Row-bind (combine) the results
ms_tab = pd.concat([pos, neg], axis=0, ignore_index=True)

# Get unique values from the 'originalID' column
ms_tab = list(ms_tab['originalID'].unique())
ms_tab.sort()

# Display the result
len(ms_tab)

# For Gene Ontology
# Run ORA with a custom GMT file
enr = gp.enrichr(
    gene_list=ms_tab, 
    gene_sets="data/EnrichR/merged_GO_Biological_Process_2021.gmt",  # Path to your custom GMT file
    organism="Human", 
    outdir="data/EnrichR/merged_GO_Biological_Process_2021_res",
    cutoff=0.1
)

results_go = enr.results

results_go_filt = results_go[(results_go['P-value'] < 0.05)]
len(results_go_filt)

results_go_filt['Genes'] = results_go_filt['Genes'].str.split(';')

results_go_filt_n = results_go_filt.explode("Genes").pivot_table(index="Term", columns="Genes", aggfunc="size", fill_value=0).reset_index()
results_go_filt_n = results_go_filt_n.set_index('Term')

results_go_filt_n.to_excel("PASNet/Input/GO/pt_fixed.xlsx")

# For KEGG
# Run ORA with a custom GMT file
enr = gp.enrichr(
    gene_list=ms_tab, 
    gene_sets="data/EnrichR/merged_KEGG_2021_Human.gmt",  # Path to your custom GMT file
    organism="Human", 
    outdir="data/EnrichR/merged_KEGG_2021_Human_res",
    cutoff=0.1
)

results_kegg = enr.results

results_kegg_filt = results_kegg[(results_kegg['P-value'] < 0.05)]
len(results_kegg_filt)

# For REACTOME
# Run ORA with a custom GMT file
enr = gp.enrichr(
    gene_list=ms_tab, 
    gene_sets="data/EnrichR/merged_Reactome_2022.gmt",  # Path to your custom GMT file
    organism="Human", 
    outdir="data/EnrichR/merged_Reactome_2022_res",
    cutoff=0.1
)

results_reactome = enr.results

results_reactome_filt = results_reactome[(results_reactome['P-value'] < 0.05)]
len(results_reactome_filt)

# For Wikipathways
# Run ORA with a custom GMT file
enr = gp.enrichr(
    gene_list=ms_tab, 
    gene_sets="data/EnrichR/merged_WikiPathway_2021_Human.gmt",  # Path to your custom GMT file
    organism="Human", 
    outdir="data/EnrichR/merged_WikiPathway_2021_res",
    cutoff=0.1
)
results_wiki = enr.results

results_wiki_filt = results_wiki[(results_wiki['P-value'] < 0.05)]
len(results_wiki_filt)

df_combined = pd.concat([results_kegg_filt, results_reactome_filt, results_wiki_filt], axis=0, ignore_index=True)

df_combined['Genes'] = df_combined['Genes'].str.split(';')

df_combined_n = df_combined.explode("Genes").pivot_table(index="Term", columns="Genes", aggfunc="size", fill_value=0).reset_index()
df_combined_n = df_combined_n.set_index('Term')

df_combined_n.to_excel("PASNet/Input/REST/pt_fixed.xlsx")

# Load activity data
activity = pd.read_csv(working_dir + "data/ICGC/activity_matrix_pvalfilt.csv")
activity['Unnamed: 0'] = activity['Unnamed: 0'].str.replace('_TF', '').str.replace('_SIG', '').str.replace('_LNC', '')
activity = activity.drop_duplicates(subset=['Unnamed: 0']).set_index('Unnamed: 0').T

# Load metadata
metadata = pd.read_csv(working_dir + "data/ICGC/metadata.csv", index_col=0)
metadata['IGHV'] = metadata['IGHV'].replace({'U': '1', 'M': '0'})

# Load pathway data
pt_rest = pd.read_excel(working_dir + "PASNet/Input/REST/pt_fixed.xlsx", index_col=0)
pt_go = pd.read_excel(working_dir + "PASNet/Input/GO/pt_fixed.xlsx", index_col=0)

# Subset activity data
activity_go = activity.loc[:, activity.columns.isin(pt_go.columns)]
activity_rest = activity.loc[:, activity.columns.isin(pt_rest.columns)]

# Merge with metadata
activity_go = activity_go.merge(metadata, left_index=True, right_index=True)
activity_rest = activity_rest.merge(metadata, left_index=True, right_index=True)

# Split into train and test sets
pretrain_go, validation_go = train_test_split(activity_go, test_size=0.2, random_state=123, stratify=activity_go['IGHV'])
pretrain_rest, validation_rest = train_test_split(activity_rest, test_size=0.2, random_state=123, stratify=activity_rest['IGHV'])

train_go, validation_grinding_go = train_test_split(pretrain_go, test_size=0.2, random_state=123, stratify=pretrain_go['IGHV'])
train_rest, validation_grinding_rest = train_test_split(pretrain_rest, test_size=0.2, random_state=123, stratify=pretrain_rest['IGHV'])

# Save to Excel
train_go.to_excel(working_dir + "PASNet/Input/GO/ICGC_training.xlsx", index=True)
validation_go.to_excel(working_dir + "PASNet/Input/GO/ICGC_validation.xlsx", index=True)
validation_grinding_go.to_excel(working_dir + "PASNet/Input/GO/ICGC_validation_grinding.xlsx", index=True)
train_rest.to_excel(working_dir + "PASNet/Input/REST/ICGC_training.xlsx", index=True)
validation_rest.to_excel(working_dir + "PASNet/Input/REST/ICGC_validation.xlsx", index=True)
validation_grinding_rest.to_excel(working_dir + "PASNet/Input/REST/ICGC_validation_grinding.xlsx", index=True)

# BCMO Data Processing
activity = pd.read_csv(working_dir + "data/ICGC/activity_matrix_pvalfilt.csv")
activity['Unnamed: 0'] = activity['Unnamed: 0'].str.replace('_TF', '').str.replace('_SIG', '').str.replace('_LNC', '')
activity = activity.drop_duplicates(subset=['Unnamed: 0']).set_index('Unnamed: 0').T

# Load BCMO metadata
metadata = pd.read_csv(working_dir + "data/BCMO/metadata.csv", index_col=0)[['IGHV']]
metadata['IGHV'] = metadata['IGHV'].replace({'U': '1', 'M': '0'})

# Subset activity data
activity_go = activity.loc[:, activity.columns.isin(pt_go.columns)]
activity_rest = activity.loc[:, activity.columns.isin(pt_rest.columns)]

# Merge with metadata
activity_go = activity_go.merge(metadata, left_index=True, right_index=True)
activity_rest = activity_rest.merge(metadata, left_index=True, right_index=True)

# Save BCMO test data
activity_go.to_excel(working_dir + "PASNet/Input/GO/BCMO_testing.xlsx", index=True)
activity_rest.to_excel(working_dir + "PASNet/Input/REST/BCMO_testing.xlsx", index=True)

