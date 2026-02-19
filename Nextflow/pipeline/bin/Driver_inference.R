#!/usr/bin/env Rscript

############################
#        LIBRARIES
############################
suppressPackageStartupMessages({
  library(dplyr)
  library(NetBID2)
  library(optparse)
})

############################
#    COMMAND LINE OPTIONS
############################
option_list <- list(
  make_option(c("-d", "--project_main_dir"), type="character", default="./",
              help="Main directory for the project [default = %default]"),
  
  make_option(c("-p", "--project_name"), type="character", default="NetBID2_Project",
              help="Project name [default = %default]"),
  
  make_option(c("--sit0"), type="character", default="M",
              help="Control group label [default = %default]"),
  
  make_option(c("--sit1"), type="character", default="U",
              help="Experiment group label [default = %default]"),
  
  make_option(c("-c", "--condition"), type="character", default="IGHV",
              help="Phenotype column name defining condition [default = %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

############################
#    ASSIGN PARAMETERS
############################
project_main_dir <- opt$project_main_dir
project_name     <- opt$project_name
sit_0            <- opt$sit0
sit_1            <- opt$sit1
Condition        <- opt$condition

driver_output <- file.path(project_main_dir, project_name, "Driver_output")

cat("========================================\n")
cat("Project:", project_name, "\n")
cat("Condition:", Condition, "\n")
cat("Comparison:", sit_1, "vs", sit_0, "\n")
cat("========================================\n")

############################
#    DIRECTORY SETUP
############################
source("bin/Adjusted_functions/NetBID.analysis.dir.create_nc.R")

analysis.par <- NetBID.analysis.dir.create_nc(
  project_main_dir = driver_output,
  project_name = project_name,
  network_dir = project_name,
  network_project_name = project_name
)

############################
# STEP 1: LOAD EXPRESSION
############################
load(file.path(project_name, "DATA/network.par.Step.exp-QC.RData"))

analysis.par$cal.eset <- network.par$net.eset

############################
# STEP 2: LOAD NETWORKS
############################
analysis.par$tf.network  <- get.SJAracne.network(network_file = analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file = analysis.par$sig.network.file)
analysis.par$nc.network  <- get.SJAracne.network(network_file = analysis.par$nc.network.file)

source("bin/Adjusted_functions/merge_TF_SIG_NC.network.R")

analysis.par$merge.network <- merge_TF_SIG_NC.network(
  TF_network  = analysis.par$tf.network,
  SIG_network = analysis.par$sig.network,
  NC_network  = analysis.par$nc.network
)

############################
# STEP 3: CALCULATE ACTIVITY
############################
ac_mat <- cal.Activity(
  target_list = analysis.par$merge.network$target_list,
  cal_mat = exprs(analysis.par$cal.eset),
  es.method = "weightedmean"
)

analysis.par$merge.ac.eset <- generate.eset(
  exp_mat = ac_mat,
  phenotype_info = pData(analysis.par$cal.eset)[colnames(ac_mat), ],
  feature_info = NULL,
  annotation_info = "activity in net-dataset"
)

############################
# DIFFERENTIAL ANALYSIS
############################
analysis.par$DE <- list()
analysis.par$DA <- list()

phe_info <- pData(analysis.par$cal.eset)

if (!(Condition %in% colnames(phe_info))) {
  stop(paste("Condition column", Condition, "not found in phenotype data"))
}

G1 <- rownames(phe_info)[phe_info[[Condition]] == sit_1]
G0 <- rownames(phe_info)[phe_info[[Condition]] == sit_0]

if (length(G1) == 0 || length(G0) == 0) {
  stop("One of the groups has zero samples. Check sit0/sit1 labels.")
}

DE_gene_bid <- getDE.BID.2G(
  eset = analysis.par$cal.eset,
  G1 = G1,
  G0 = G0,
  G1_name = sit_1,
  G0_name = sit_0
)

DA_driver_bid <- getDE.BID.2G(
  eset = analysis.par$merge.ac.eset,
  G1 = G1,
  G0 = G0,
  G1_name = sit_1,
  G0_name = sit_0
)

comparison_name <- paste0(sit_1, ".Vs.", sit_0)

analysis.par$DE[[comparison_name]] <- DE_gene_bid
analysis.par$DA[[comparison_name]] <- DA_driver_bid

NetBID.saveRData(analysis.par = analysis.par, step = "act-DA")

############################
# STEP 4: MASTER TABLE
############################
NetBID.loadRData(analysis.par = analysis.par, step = "act-DA")

db.preload(use_level = "gene", use_spe = "human", update = FALSE)

all_comp <- names(analysis.par$DE)

use_genes <- unique(c(
  analysis.par$merge.network$network_dat$source.symbol,
  analysis.par$merge.network$network_dat$target.symbol
))

analysis.par$transfer_tab <- get_IDtransfer2symbol2type(
  from_type = "hgnc_symbol",
  use_genes = use_genes
)

source("bin/Adjusted_functions/generate.masterTable_nc.R")

analysis.par$final_ms_tab <- generate.masterTable_lnc(
  use_comp = all_comp,
  DE = analysis.par$DE,
  DA = analysis.par$DA,
  target_list = analysis.par$merge.network$target_list,
  tf_sigs = tf_sigs,
  z_col = "Z-statistics",
  display_col = c("logFC", "adj.P.Val"),
  main_id_type = "hgnc_symbol"
)

############################
# EXPORT RESULTS
############################
out_file <- file.path(
  analysis.par$out.dir.DATA,
  paste0(analysis.par$project.name, "_ms_tab.xlsx")
)

out2excel(analysis.par$final_ms_tab, out.xlsx = out_file)

NetBID.saveRData(analysis.par = analysis.par, step = "ms-tab")

cat("Analysis completed successfully.\n")
cat("Output saved to:", out_file, "\n")
