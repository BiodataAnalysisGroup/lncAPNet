#!/usr/bin/python

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