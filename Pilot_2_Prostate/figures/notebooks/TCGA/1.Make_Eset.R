
library(tidyverse)
library(Biobase)
library(umap)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(openxlsx)

rm(list = ls())
gc()

# load data

metadata <- read.xlsx("../initial_data/metadata.xlsx")
count_matrix <- read.xlsx("../initial_data/count_matrix.xlsx")

# merge metadata

rownames(metadata) <- metadata$submitter_id
metadata_no_normal <- filter(metadata, Status != "Normal" & Status != "Paracancerous")
metadata_no_normal <- metadata_no_normal %>%
  mutate(Sample_ID = str_remove(Sample_ID, "-01A.*"))
rownames(metadata_no_normal) <- metadata_no_normal$Sample_ID

metadata_no_normal <- merge(metadata_no_normal, metadata, by = 0)
metadata_no_normal <- metadata_no_normal[!is.na(metadata_no_normal$ajcc_pathologic_n), ]
rownames(metadata_no_normal) <- metadata_no_normal$Row.names

## fix name in count matrix

list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol', 'gene_biotype'), 
              filters ="ensembl_gene_id",
              value = count_matrix$X1,
              mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
)

list <- list[list$hgnc_symbol != "", ]
list <- list[!duplicated(list$ensembl_gene_id),]
list <- list[!duplicated(list$hgnc_symbol),]
count_matrix <- count_matrix[count_matrix$X1 %in% list$ensembl_gene_id,]

rownames(count_matrix) <- count_matrix$X1
rownames(list) <- list$ensembl_gene_id
count_matrix <- merge(count_matrix, list, by = 0)
rownames(count_matrix) <- count_matrix$hgnc_symbol
count_matrix <- subset(count_matrix, select = -c(Row.names, X1, ensembl_gene_id, hgnc_symbol, gene_biotype))
rownames(list) <- list$hgnc_symbol

## count_matrix pre-processing

count_matrix <- as.data.frame(t(count_matrix))
count_matrix$samples <- rownames(count_matrix)

count_matrix <- count_matrix %>%
  mutate(samples = str_remove(samples, "-01A.*"))

count_matrix <- count_matrix[order(count_matrix$samples),]
rownames(count_matrix) <- count_matrix$samples
count_matrix <- subset(count_matrix, select = -c(samples))
count_matrix <- as.data.frame(t(count_matrix))

count_matrix_no_normal <- count_matrix[,colnames(count_matrix) %in% metadata_no_normal$Sample_ID]

# Transpose count matrix: rows = samples, columns = genes
data_pca <- t(count_matrix_no_normal)

# Perform PCA (center and scale the data)
pca_result <- prcomp(data_pca, center = TRUE, scale. = FALSE)

# Create PCA dataframe with PC1 and PC2
pca_df <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(data_pca)

# Add Status from metadata
metadata_no_normal$Sample <- rownames(metadata_no_normal)
pca_df <- left_join(pca_df, metadata_no_normal, by = "Sample")

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = ajcc_pathologic_n)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Samples", color = "Status")



# Transpose count matrix: samples as rows, genes as columns
data_umap <- t(count_matrix_no_normal)

# Run UMAP
umap_result <- umap(data_umap)

# Create UMAP dataframe
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample_ID <- rownames(data_umap)

# Add Status from metadata (assuming rownames(metadata) == colnames(count_matrix))
metadata_no_normal$Sample_ID <- rownames(metadata_no_normal)
umap_df <- merge(umap_df, metadata_no_normal, by = 0)

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = ajcc_pathologic_n)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP of Samples", color = "Status")


all(colnames(count_matrix_no_normal)==rownames(metadata_no_normal))
all(rownames(count_matrix_no_normal)==rownames(list))

count_matrix <- as.matrix(count_matrix_no_normal)

ExpressionSet <- ExpressionSet(assayData = count_matrix, phenoData = new("AnnotatedDataFrame", data = metadata_no_normal), featureData = new("AnnotatedDataFrame", data = list))
saveRDS(ExpressionSet, file = "../data/net_eset.rds")
