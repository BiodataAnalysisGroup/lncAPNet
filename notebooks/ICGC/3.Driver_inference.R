#!/usr/bin/R

# Load necessary libraries
library(dplyr)
library(NetBID2)
library(optparse)

rm(list = ls())
gc()

setwd("C:/Users/vasileioubill95/Desktop/Projects/lncAPNet_total/")

# Assign input values to script variables

project_main_dir <- './'  # user defined main directory for the project
project_name <- "ICGC_NetBID2" # project name for the project folders under main directory
sit_0 <- "M"
sit_1 <- "U"
Condition <- "IGHV"

driver_output <- sprintf('%s/Driver_output', project_name)

source("scripts/Adjusted_functions/NetBID.analysis.dir.create_nc.R")
# This list object (analysis.par) is an ESSENTIAL variable in driver estimation pipeline
analysis.par  <- NetBID.analysis.dir.create_nc(project_main_dir=driver_output, project_name=project_name,
                                                network_dir=project_name, network_project_name=project_name)

############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############

# If use the same expression dataset as in the network construction, just reload it directly
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData', project_name)) # RData saved after QC in the network construction step
analysis.par$cal.eset <- network.par$net.eset

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$nc.network <- get.SJAracne.network(network_file=analysis.par$nc.network.file)

# # Create QC report for the network
# Merge network first
source("scripts/Adjusted_functions/merge_TF_SIG_NC.network.R")
analysis.par$merge.network <- merge_TF_SIG_NC.network(TF_network=analysis.par$tf.network,
                                                       SIG_network=analysis.par$sig.network,
                                                       NC_network = analysis.par$nc.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info[[Condition]]==sit_1)] # Experiment group
G0  <- rownames(phe_info)[which(phe_info[[Condition]]==sit_0)] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name=sit_1,G0_name=sit_0)
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name=sit_1,G0_name=sit_0)
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[['U.Vs.M']] <- DE_gene_bid
analysis.par$DA[['U.Vs.M']] <- DA_driver_bid

# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

############### Step 4: Generate a master table for drivers (ms-tab) ###############

# Reload analysis.par RData from Step 3
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'hgnc_symbol',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Creat the final master table
source("scripts/Adjusted_functions/generate.masterTable_nc.R")
analysis.par$final_ms_tab <- generate.masterTable_lnc(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                      target_list=analysis.par$merge.network$target_list,
                                                      tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','adj.P.Val'),
                                                      main_id_type='hgnc_symbol')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)

# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file)

# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')
