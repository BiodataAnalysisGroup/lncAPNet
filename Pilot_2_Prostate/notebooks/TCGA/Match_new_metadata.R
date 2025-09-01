
library(tidyverse)
library(Biobase)
library(openxlsx)

eset <- readRDS("../data/net_eset.rds")
pheno <- pData(eset)
pheno_new <- read.xlsx("../initial_data/PRAD_Clin.xlsx")
rownames(pheno_new) <- pheno_new$submitter_id

df_clean <- pheno %>%
  mutate(name = str_remove(Sample_ID, "-01A.*"))

rownames(df_clean) <- df_clean$name

df <- merge(df_clean, pheno_new, by = 0)
rownames(df) <- df$Sample_ID

expression <- as.data.frame(exprs(eset))
expression <- expression[,colnames(expression) %in% df$Sample_ID]
expression <- as.data.frame(t(expression))
expression <- expression[order(rownames(expression)),]
expression <- as.data.frame(t(expression))
expression <- as.matrix(expression)

features <- fData(eset)

all(colnames(expression)==rownames(df))
all(rownames(expression)==rownames(features))


ExpressionSet <- ExpressionSet(assayData = expression, phenoData = new("AnnotatedDataFrame", data = df), featureData = new("AnnotatedDataFrame", data = features))
