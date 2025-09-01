#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(Biobase)
library(optparse)
library(biomaRt)

rm(list = ls())  # Clear workspace
gc()  # Garbage collection to free up memory

setwd("C:/Users/vasileioubill95/Desktop/Projects/lncAPNet_total/")

# Input data
counts_table <- read.csv("data/BCMO/count_matrix.csv", row.names = 1)
metadata <- read.csv("data/BCMO/metadata.csv", row.names = 1)

# Ensure columns in counts_table match rownames in metadata
counts_table <- counts_table[, colnames(counts_table) %in% rownames(metadata)]

list <- getBM(attributes=c('hgnc_symbol', 'gene_biotype'), 
              filters ='hgnc_symbol',
              value = rownames(counts_table),
              mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
)

table(list$gene_biotype)
list_mrnas <- filter(list, gene_biotype == "protein_coding")
list_nc <- filter(list, gene_biotype != "protein_coding")

write.table(list_nc$hgnc_symbol, file = "data/BCMO/nc_list.txt", row.names = F, quote = F)

list_mrnas_nc <- rbind(list_mrnas, list_nc)

counts_table <- counts_table[rownames(counts_table) %in% list_mrnas_nc$hgnc_symbol,]

# Create ExpressionSet
eset <- ExpressionSet(assayData = as.matrix(counts_table), phenoData = AnnotatedDataFrame(metadata))

# Save the ExpressionSet object
saveRDS(eset, "data/BCMO/net_eset.rds")

