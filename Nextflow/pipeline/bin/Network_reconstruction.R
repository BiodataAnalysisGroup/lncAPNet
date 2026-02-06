#!/usr/bin/env Rscript

###############################################
### Step 0: Parse command-line arguments first
###############################################

Sys.setenv(
  XDG_CACHE_HOME = file.path(tempdir(), "r_cache")
)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--eset"), type="character", help="Path to input RDS eset file [REQUIRED]"),
  make_option(c("--project_dir"), type="character", default="./",
              help="Main project directory [default=./]"),
  make_option(c("--project_name"), type="character", default="NetBID2_Project",
              help="Project name [default=NetBID2_Project]"),
  make_option(c("--iqr"), type="double", default=0.5,
              help="IQR threshold for whole genes [default=0.5]"),
  make_option(c("--iqr_loose"), type="double", default=0.1,
              help="IQR loose threshold for hub genes [default=0.1]"),
  make_option(c("--log"), type="logical", default=TRUE,
              help="Print log messages for normalization [default=TRUE]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# If --eset is missing, print help and exit
if (is.null(opt$eset)) {
  print_help(opt_parser)
  stop("Error: --eset argument is required.", call.=FALSE)
}

###############################################
### Step 1: Load heavy libraries AFTER parsing args
###############################################
suppressPackageStartupMessages({
  library(dplyr)
  library(NetBID2)
  library(biomaRt)
})

###############################################
### Step 2: Load eset
###############################################
message("Loading eset: ", opt$eset)
net_eset <- readRDS(opt$eset)

###############################################
### Step 3: Set project directories
###############################################
project_main_dir <- opt$project_dir
project_name     <- opt$project_name

network.par <- NetBID.network.dir.create(
  project_main_dir = project_main_dir,
  project_name     = project_name
)

network.par$net.eset <- net_eset

###############################################
### Step 4: Normalization (conditional log2)
###############################################
if (opt$log) {
  message("Applying log2 normalization...")
  mat <- log2(exprs(network.par$net.eset) + 1)
} else {
  mat <- exprs(network.par$net.eset)
}

message("Filtering low-expression genes...")
choose1 <- apply(mat <= quantile(mat, probs = 0.05), 1, sum) <= ncol(mat) * 0.90
mat <- mat[choose1, ]

net_eset <- generate.eset(
  exp_mat        = mat,
  phenotype_info = pData(network.par$net.eset)[colnames(mat), ],
  feature_info   = fData(network.par$net.eset)[rownames(mat), ],
  annotation_info= annotation(network.par$net.eset)
)

network.par$net.eset <- net_eset

message("Saving QC data...")
NetBID.saveRData(network.par = network.par, step = "exp-QC")

###############################################
### Step 5: Database + TF/SIG selection
###############################################
message("Loading TF/SIG database...")
db.preload(use_level='gene', use_spe='human', update=FALSE)

use_gene_type <- 'hgnc_symbol'
use_genes <- rownames(fData(network.par$net.eset))
use_list <- get.TF_SIG.list(use_genes, use_gene_type = use_gene_type)

###############################################
### Step 6: Prepare SJAracne
###############################################
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe)
prj.name <- network.par$project.name

message("Preparing SJAracne...")

SJAracne.prepare(
  eset              = network.par$net.eset,
  use.samples       = use.samples,
  TF_list           = use_list$tf,
  SIG_list          = use_list$sig,
  IQR.thre          = opt$iqr,
  IQR.loose_thre    = opt$iqr_loose,
  SJAR.project_name = prj.name,
  SJAR.main_dir     = network.par$out.dir.SJAR
)

###############################################
### Step 7: Retrieve non-protein-coding genes (lncRNAs) using biomaRt
###############################################
message("Retrieving non-protein-coding genes from biomaRt...")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters    = "hgnc_symbol",
  values     = rownames(fData(network.par$net.eset)),
  mart       = ensembl
)

lncRNAs <- gene_info %>%
  filter(gene_biotype != "protein_coding") %>%
  pull(hgnc_symbol)

message("Retrieved ", length(lncRNAs), " non-protein-coding genes (lncRNAs).")

###############################################
### Step 8: Filtering lncRNA list
###############################################
message("Filtering lncRNA entries ...")
input_exp_file <- file.path(project_main_dir, project_name,
                            "SJAR", project_name, "input.exp")
exp_table <- read.table(input_exp_file, header = TRUE)

lnc_filtered <- lncRNAs[lncRNAs %in% exp_table$isoformId]

out_lnc_file <- file.path(project_main_dir, project_name,
                          "SJAR", project_name, "lnc_list.txt")

write.table(lnc_filtered, out_lnc_file,
            row.names = FALSE, quote = FALSE, col.names = FALSE)

message("Filtered lncRNA list saved to: ", out_lnc_file)

###############################################
### End
###############################################
message("Done.")
