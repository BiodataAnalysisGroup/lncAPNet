library(tidyverse)
library(openxlsx)
library(biomaRt)

rm(list = ls())
gc()

enr <- read.xlsx("C:/Users/vasileioubill95/Desktop/Projects/lncAPNet_Prostate/data/EnrichR/both/merged_GO_Biological_Process_2021.xlsx")
enr <- subset(enr, select = -c(Description))

ms_tab <- read.xlsx("C:/Users/vasileioubill95/Desktop/Projects/lncAPNet_Prostate/Prostate_NetBID2/Driver_output/Prostate_NetBID2/DATA/ms_tab_gene_annotation_fixed.xlsx")


changes <- 0  # counter for replacements

for (col in names(enr)[-1]) {
  print(col)
  enr_updated[[col]] <- sapply(enr[[col]], function(gene) {
    match_idx <- which(list_features$ensembl_gene_id == gene)
    if (length(match_idx) > 0) {
      changes <<- changes + 1  # increment global count of changes
      return(list_features$hgnc_symbol[match_idx[1]])
    } else {
      return(gene)
    }
  })
}

cat("Total replacements made:", changes, "\n")

# Optionally, show rows/columns where changes happened
diff_positions <- which(enr != enr_updated, arr.ind = TRUE)
if (nrow(diff_positions) > 0) {
  print("Positions (row, col) where replacements happened:")
  print(diff_positions)
} else {
  cat("No changes detected.\n")
}

enr_updated <- read.xlsx("C:/Users/vasileioubill95/Desktop/test.xlsx")




# Ensure the dataframes have the same structure
if (!all(names(enr) == names(enr_updated)) || nrow(enr) != nrow(enr_updated)) {
  stop("Dataframes must have the same dimensions and column names.")
}

# Create a function to compare cell values
compare_dfs <- function(df1, df2) {
  differences <- data.frame(
    row = integer(),
    column = character(),
    old_value = character(),
    new_value = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(df1))) {
    for (j in seq_along(df1)) {
      val1 <- as.character(df1[i, j])
      val2 <- as.character(df2[i, j])
      
      if (!identical(val1, val2)) {
        differences <- rbind(
          differences,
          data.frame(
            row = i,
            column = names(df1)[j],
            old_value = val1,
            new_value = val2,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  return(differences)
}

# Run the comparison
differences_list <- compare_dfs(enr, enr_updated)

# View differences
print(differences_list)

