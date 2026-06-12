library(tidyverse)
library(openxlsx)
library(scMINER)
library(anndata)
library(icesTAF)
library(Seurat)
library(reticulate)
library(SingleCellExperiment)
library(biomaRt)

rm(list = ls())
gc()

adata <- read_h5ad("../../../data/TCGA_Breast/subset_epithelial_data_total.h5ad")
adata_seurat <- CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)

saveRDS(adata_seurat, file = "../../../data/TCGA_Breast/brca_epithelial_subset_seurat_one_network.rds")

rm(list = ls())
gc()

## Define Directory

setwd("./")


## Load Seurat Object

SeuratObject <- readRDS("../../../data/TCGA_Breast/brca_epithelial_subset_seurat_one_network.rds")

scminer_dir <- createProjectSpace(project_dir = "./", project_name = "scRNAseq_Breast_Atlas_one_network")

## Convert Seurat Object into SparseEset

meta.data <- SeuratObject@meta.data
meta.data$celltype_minor <- droplevels(meta.data$subtype_pure)
table(meta.data$subtype_pure)

feature.data<-data.frame(rownames(SeuratObject@assays$RNA))
colnames(feature.data)<-"geneSymbol"
rownames(feature.data)<-feature.data$geneSymbol

expression <- as.data.frame(SeuratObject@assays$RNA$counts)
expression <- expression[,colnames(expression) %in% rownames(meta.data)]


eset <- createSparseEset(input_matrix =expression, 
                         cellData = meta.data, 
                         featureData =  feature.data, 
                         projectID = "scRNAseq_Breast_Atlas_BasalvsLum_one_network", 
                         addMetaData = TRUE)

saveRDS(eset, file = "../../../data/TCGA_Breast/brca_epithelial_subset_seurat_eset_one_network.rds")

eset <- readRDS("../../../data/TCGA_Breast/brca_epithelial_subset_seurat_eset_one_network.rds")

## Create SJARACNe input Objects by taking into consideration:
##    a) transcription factors (TFs) as gene hubs
##    b) Signalling genes (SIGs) as gene hubs

## Columns with any illegal characters can not be used for groupping
generateSJARACNeInput(input_eset = eset, group_name = "normal_cell_call", 
                      sjaracne_dir = "./scRNAseq_Breast_Atlas_one_network/SJARACNe", 
                      species_type = "hg", 
                      driver_type = "TF_SIG",
                      downSample_N = 2000)

genes <- rownames(as.data.frame(exprs(eset)))

mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl")

list_type <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"), 
  filters = "hgnc_symbol",                       
  values = genes,                           
  mart = mart
)


print(list_type)

lnc_list <- list_type[list_type$gene_biotype != "protein_coding",]

case <- read.table("scRNAseq_Breast_Atlas_one_network/SJARACNe/cancer/cancer.22712_2000.exp.txt", header = T)

lnc_list <- lnc_list[lnc_list$hgnc_symbol %in% case$isoformId,]
lnc_list <- lnc_list$hgnc_symbol

dir.create("scRNAseq_Breast_Atlas_one_network/SJARACNe/cancer/LNC", recursive = TRUE)
write.table(lnc_list, file = "scRNAseq_Breast_Atlas_one_network/SJARACNe/cancer/LNC/lnc_list.txt", row.names = F, quote = F, col.names = F)

# sed 's/[[:space:]]*$//' sig.txt > sig_cleaned.txt before SJARACNe
# sed 's/[[:space:]]*$//' tf.txt > tf_cleaned.txt before SJARACNe
# sed 's/[[:space:]]*$//' lnc_list.txt > lnc_cleaned.txt before SJARACNe

## RUN SJARACNe

# run bash basic Commands with SJARACNe algorithms 


## Calculate Activity

eset <- readRDS("../../../data/TCGA_Breast/brca_epithelial_subset_seurat_eset.rds")

source("../Adjusted_functions/getActivity_inBatch_lnc.R")

activity.eset <- getActivity_inBatch_lnc(input_eset = eset,
                                     sjaracne_dir = "./scRNAseq_Breast_Atlas_one_network/SJARACNe",
                                     group_name = "normal_cell_call", 
                                     driver_type = "TF_SIG_LNC", 
                                     activity_method = "mean", 
                                     do.z_normalization = F)

saveRDS(activity.eset, file = "../../../data/TCGA_Breast/brca_epithelial_subset_seurat_eset_activity.rds")

activity.eset <- readRDS("../../../data/TCGA_Breast/brca_epithelial_subset_seurat_eset_activity.rds")
activity <- as.data.frame(exprs(activity.eset))

metadata <- pData(activity.eset)

write.csv(activity, file = "../../../data/TCGA_Breast/brca_epithelial_activity_matrix.csv", row.names = T, quote = F)

write.csv(metadata, file = "../../../data/TCGA_Breast/brca_epithelial_metadata.csv", row.names = T, quote = F)
