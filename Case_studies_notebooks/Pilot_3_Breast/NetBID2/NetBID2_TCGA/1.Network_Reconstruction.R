#!/usr/bin/R

# Load necessary libraries
library(dplyr)
library(NetBID2)
library(optparse)
library(org.Hs.eg.db)

# Clear workspace and perform garbage collection
rm(list = ls())  
gc()

setwd("C:/Users/vvasileiou/Desktop/test_lncAPNet/Pilot_3_TCGA_Breast/notebooks/TCGA_Breast/")

# Load RDS eset file
net_eset <- readRDS("../../../data/TCGA_BRCA_raw_ExpressionSet.rds")

# Define main working directory and project name
project_main_dir <- './'  # user defined main directory for the project
project_name <- "Breast_NetBID2" # project name for the project folders under main directory

# This list object (network.par) is an ESSENTIAL variable in the network construction pipeline
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir, project_name=project_name)

network.par$net.eset <- net_eset

# Get the expression matrix from ExpressionSet object
mat <- exprs(network.par$net.eset)

# Apply log2 transformation if --norm is TRUE
mat <- log2(mat + 1)

# Perform limma quantile normalization
# Filter out genes with very low expression values (bottom 5%) in most samples (more than 90%)
choose1 <- apply(mat <= quantile(mat, probs = 0.05), 1, sum) <= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]

# Update eset with normalized expression matrix
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))

network.par$net.eset <- net_eset

NetBID.saveRData(network.par = network.par,step='exp-QC')

# Load the database
db.preload(use_level='gene', use_spe='human', update=FALSE)

# Convert gene ID into the corresponding TF/SIG list
use_gene_type <- 'hgnc_symbol' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes, use_gene_type=use_gene_type)

# Select samples for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # using all samples, can modify if needed
prj.name <- network.par$project.name # if using different samples, change the project name

# Prepare and execute the SJAracne algorithm with user-defined IQR thresholds
SJAracne.prepare(
  eset=network.par$net.eset,
  use.samples=use.samples,
  TF_list=use_list$tf,
  SIG_list=use_list$sig,
  IQR.thre = 0.5,        # IQR threshold for whole genes
  IQR.loose_thre = 0.1,  # IQR threshold for hub genes
  SJAR.project_name=prj.name,
  SJAR.main_dir=network.par$out.dir.SJAR
)

table <- read.table("Breast_NetBID2/SJAR/Breast_NetBID2/input.exp", header = T)

feature <- fData(net_eset)
lnc <- feature[feature$gene_biotype !="protein_coding",]
lnc <- lnc$hgnc_symbol

lnc_list_filt <- lnc[lnc %in% table$isoformId]

# Write the concatenated content to a new file
write.table(lnc_list_filt, "Breast_NetBID2/SJAR/Breast_NetBID2/lnc_list.txt", row.names = F, quote = F, col.names = F)

# sed 's/[[:space:]]*$//' sig.txt > sig_cleaned.txt before SJARACNe
# sed 's/[[:space:]]*$//' tf.txt > tf_cleaned.txt before SJARACNe
# sed 's/[[:space:]]*$//' lnc_list.txt > lnc_cleaned.txt before SJARACNe