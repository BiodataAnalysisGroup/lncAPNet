#!/usr/bin/env Rscript

############################
#        LIBRARIES
############################
suppressPackageStartupMessages({
  library(dplyr)
  library(NetBID2)
  library(optparse)
})

## Adjust function of NetBID.analysis.dir.create_nc

NetBID.analysis.dir.create_nc <- function (project_main_dir = NULL, project_name = NULL, network_dir = NULL, 
          network_project_name = NULL, tf.network.file = NULL, sig.network.file = NULL) 
{
  analysis.par <- list()
  analysis.par$main.dir <- project_main_dir
  analysis.par$project.name <- project_name
  analysis.par$out.dir <- sprintf("%s/%s/", analysis.par$main.dir, 
                                  analysis.par$project.name)
  analysis.par$tf.network.file <- ""
  analysis.par$sig.network.file <- ""
  analysis.par$nc.network.file <- ""
  # TF
  if (is.null(tf.network.file) == FALSE) {
    analysis.par$tf.network.file <- tf.network.file
  }
  else {
    tf_net1 <- sprintf("%s/SJAR/%s/output_tf_sjaracne_%s_out_.final/consensus_network_ncol_.txt", 
                       network_dir, network_project_name, network_project_name)
    tf_net2 <- sprintf("%s/SJAR/SJARACNE_%s_TF/consensus_network_ncol_.txt", 
                       network_dir, network_project_name)
    if (file.exists(tf_net2)) 
      analysis.par$tf.network.file <- tf_net2
    else analysis.par$tf.network.file <- tf_net1
  }
  # SIG
  if (is.null(tf.network.file) == FALSE) {
    analysis.par$sig.network.file <- sig.network.file
  }
  else {
    sig_net1 <- sprintf("%s/SJAR/%s/output_sig_sjaracne_%s_out_.final/consensus_network_ncol_.txt", 
                        network_dir, network_project_name, network_project_name)
    sig_net2 <- sprintf("%s/SJAR/SJARACNE_%s_SIG/consensus_network_ncol_.txt", 
                        network_dir, network_project_name)
    if (file.exists(sig_net2)) 
      analysis.par$sig.network.file <- sig_net2
    else analysis.par$sig.network.file <- sig_net1
  }
  # NC
  if (is.null(tf.network.file) == FALSE) {
    analysis.par$nc.network.file <- nc.network.file
  }
  else {
    nc_net1 <- sprintf("%s/SJAR/%s/output_nc_sjaracne_%s_out_.final/consensus_network_ncol_.txt", 
                        network_dir, network_project_name, network_project_name)
    nc_net2 <- sprintf("%s/SJAR/SJARACNE_%s_NC/consensus_network_ncol_.txt", 
                        network_dir, network_project_name)
    if (file.exists(sig_net2)) 
      analysis.par$nc.network.file <- nc_net2
    else analysis.par$nc.network.file <- nc_net1
  }
  #TF
  if (file.exists(analysis.par$tf.network.file)) {
    message(sprintf("TF network file found in %s", analysis.par$tf.network.file))
  }
  else {
    message(sprintf("TF network file not found in %s, please check and re-try !", 
                    analysis.par$tf.network.file))
    return(FALSE)
  }
  #SIG
  if (file.exists(analysis.par$sig.network.file)) {
    message(sprintf("SIG network file found in %s", analysis.par$sig.network.file))
  }
  else {
    message(sprintf("SIG network file not found in %s, please check and re-try ", 
                    analysis.par$sig.network.file))
    return(FALSE)
  }
  #NC
  if (file.exists(analysis.par$nc.network.file)) {
    message(sprintf("NC network file found in %s", analysis.par$nc.network.file))
  }
  else {
    message(sprintf("NC network file not found in %s, please check and re-try !", 
                    analysis.par$nc.network.file))
    return(FALSE)
  }
  if (!dir.exists(analysis.par$out.dir)) {
    dir.create(analysis.par$out.dir, recursive = TRUE)
  }
  analysis.par$out.dir.QC <- paste0(analysis.par$out.dir, 
                                    "/QC/")
  if (!dir.exists(analysis.par$out.dir.QC)) {
    dir.create(analysis.par$out.dir.QC, recursive = TRUE)
  }
  analysis.par$out.dir.DATA <- paste0(analysis.par$out.dir, 
                                      "/DATA/")
  if (!dir.exists(analysis.par$out.dir.DATA)) {
    dir.create(analysis.par$out.dir.DATA, recursive = TRUE)
  }
  analysis.par$out.dir.PLOT <- paste0(analysis.par$out.dir, 
                                      "/PLOT/")
  if (!dir.exists(analysis.par$out.dir.PLOT)) {
    dir.create(analysis.par$out.dir.PLOT, recursive = TRUE)
  }
  message(sprintf("Analysis space created, please check %s", 
                  analysis.par$out.dir))
  return(analysis.par)
}

## Adjust function of merge_TF_SIG_NC.network

merge_TF_SIG_NC.network <- function (TF_network = NULL, SIG_network = NULL, 
                                      NC_network = NULL) 
{
  s_TF <- names(TF_network$target_list)
  s_SIG <- names(SIG_network$target_list)
  s_NC <- names(NC_network$target_list)
  funcType <- c(rep("TF", base::length(s_TF)), rep("SIG",base::length(s_SIG)),
                rep("NC", base::length(s_NC)))
  rn <- c(s_TF, s_SIG, s_NC)
  rn_label <- base::paste(rn, funcType, sep = "_")
  target_list_combine <- c(TF_network$target_list, SIG_network$target_list, 
                           NC_network$target_list)
  names(target_list_combine) <- rn_label
  #TF
  n_TF <- TF_network$network_dat
  if (nrow(n_TF) > 0) 
    n_TF$source <- base::paste(n_TF$source, "TF", sep = "_")
  #SIG
  n_SIG <- SIG_network$network_dat
  if (nrow(n_SIG) > 0) 
    n_SIG$source <- base::paste(n_SIG$source, "SIG", sep = "_")
  #NC
  n_NC <- NC_network$network_dat
  if (nrow(n_NC) > 0) 
    n_NC$source <- base::paste(n_NC$source, "NC", sep = "_")
  net_dat <- base::rbind(n_TF, n_SIG, n_NC)
  igraph_obj <- graph_from_data_frame(net_dat[, c("source", 
                                                  "target")], directed = TRUE)
  if ("MI" %in% colnames(net_dat)) 
    igraph_obj <- set_edge_attr(igraph_obj, "weight", index = E(igraph_obj), 
                                value = net_dat[, "MI"])
  if ("spearman" %in% colnames(net_dat)) 
    igraph_obj <- set_edge_attr(igraph_obj, "sign", index = E(igraph_obj), 
                                value = sign(net_dat[, "spearman"]))
  return(list(network_dat = net_dat, target_list = target_list_combine, 
              igraph_obj = igraph_obj))
}

## Adjust function of generate.masterTable_lnc

generate.masterTable_lnc <- function (use_comp = NULL, DE = NULL, DA = NULL, target_list = NULL, 
          main_id_type = NULL, transfer_tab = NULL, tf_sigs = NULL, 
          z_col = "Z-statistics", display_col = c("logFC", "P.Value"), 
          column_order_strategy = "type") 
{
  ori_rn <- rownames(DA[[1]])
  w1 <- grep("(.*)_TF", ori_rn)
  w2 <- grep("(.*)_SIG", ori_rn)
  w3 <- grep("(.*)_NC", ori_rn)
  funcType <- rep(NA, length.out = base::length(ori_rn))
  rn <- funcType
  funcType[w1] <- "TF"
  funcType[w2] <- "SIG"
  funcType[w3] <- "NC"
  rn[w1] <- gsub("(.*)_TF", "\\1", ori_rn[w1])
  rn[w2] <- gsub("(.*)_SIG", "\\1", ori_rn[w2])
  rn[w3] <- gsub("(.*)_NC", "\\1", ori_rn[w3])
  rn_label <- ori_rn
  use_size <- unlist(lapply(target_list[rownames(DA[[1]])], 
                            nrow))
  current_id <- names(tf_sigs$tf)[-1]
  use_info <- base::unique(base::rbind(tf_sigs$tf$info, tf_sigs$sig$info))
  if (main_id_type %in% current_id) {
    use_info <- use_info[which(use_info[, main_id_type] %in% 
                                 rn), ]
  }
  else {
    if (is.null(transfer_tab) == TRUE) {
      transfer_tab <- get_IDtransfer(from_type = main_id_type, 
                                     to_type = current_id[1], use_genes = rn, ignore_version = TRUE)
      use_info <- base::merge(use_info, transfer_tab, 
                              by.x = current_id[1], by.y = current_id[1])
    }
    else {
      transfer_tab <- transfer_tab[which(transfer_tab[, 
                                                      main_id_type] %in% rn), ]
      uid <- base::intersect(colnames(transfer_tab), colnames(use_info))
      if (base::length(uid) == 0) {
        message("No ID type in the transfer_tab could match ID type in tf_sigs, please check and re-try!")
        return(FALSE)
      }
      uid <- uid[1]
      use_info <- base::merge(use_info, transfer_tab, 
                              by.x = uid, by.y = uid)
    }
    use_info <- use_info[which(use_info[, main_id_type] %in% 
                                 rn), ]
  }
  use_info <- base::unique(use_info)
  if (nrow(use_info) == 0) {
    message("ID issue error, please check main_id_type setting!")
    return(FALSE)
  }
  tmp1 <- stats::aggregate(use_info, list(use_info[, main_id_type]), 
                           function(x) {
                             x1 <- x[which(x != "")]
                             x1 <- x1[which(is.na(x1) == FALSE)]
                             base::paste(sort(base::unique(x1)), collapse = ";")
                           })
  tmp1 <- tmp1[, -1]
  rownames(tmp1) <- tmp1[, main_id_type]
  geneSymbol <- tmp1[rn, "external_gene_name"]
  if ("external_transcript_name" %in% colnames(tmp1)) {
    gene_label <- base::paste(tmp1[rn, "external_transcript_name"], 
                              funcType, sep = "_")
  }
  else {
    gene_label <- base::paste(tmp1[rn, "external_gene_name"], 
                              funcType, sep = "_")
  }
  label_info <- data.frame(originalID_label = rn_label, originalID = rn, 
                           gene_label = gene_label, geneSymbol = geneSymbol, funcType = funcType, 
                           Size = use_size, stringsAsFactors = FALSE)
  w1 <- which(is.na(geneSymbol) == TRUE)
  label_info[w1, "geneSymbol"] <- label_info[w1, "originalID"]
  label_info[w1, "gene_label"] <- label_info[w1, "originalID_label"]
  add_info <- tmp1[rn, ]
  combine_info <- lapply(use_comp, function(x) {
    DA[[x]] <- DA[[x]][rn_label, , drop = F]
    DE[[x]] <- as.data.frame(DE[[x]])[rn, ]
    avg_col <- colnames(DA[[x]])[grep("^Ave", colnames(DA[[x]]))]
    uc <- c(z_col, avg_col, base::setdiff(display_col, c(z_col, 
                                                         avg_col)))
    uc <- base::intersect(uc, colnames(DA[[x]]))
    DA_info <- DA[[x]][rn_label, uc, drop = F]
    avg_col <- colnames(DE[[x]])[grep("^Ave", colnames(DE[[x]]))]
    uc <- c(z_col, avg_col, base::setdiff(display_col, c(z_col, 
                                                         avg_col)))
    uc <- base::intersect(uc, colnames(DE[[x]]))
    DE_info <- as.data.frame(DE[[x]])[rn, uc, drop = F]
    colnames(DA_info) <- paste0(colnames(DA_info), ".", 
                                x, "_DA")
    colnames(DE_info) <- paste0(colnames(DE_info), ".", 
                                x, "_DE")
    colnames(DA_info)[1] <- paste0("Z.", x, "_DA")
    colnames(DE_info)[1] <- paste0("Z.", x, "_DE")
    out <- base::cbind(DA_info, DE_info, stringsAsFactors = FALSE)
    rownames(out) <- rn_label
    out
  })
  combine_info_DA <- do.call(base::cbind, lapply(combine_info, 
                                                 function(x) x[rn_label, grep("_DA$", colnames(x)), drop = T]))
  combine_info_DE <- do.call(base::cbind, lapply(combine_info, 
                                                 function(x) x[rn_label, grep("_DE$", colnames(x)), drop = T]))
  if (column_order_strategy == "type" & length(use_comp) > 
      1) {
    col_ord <- c("Z", "AveExpr", display_col)
    tmp1 <- lapply(col_ord, function(x) {
      x1 <- grep(sprintf("^%s\\.", x), colnames(combine_info_DA))
      if (length(x1) > 0) 
        combine_info_DA[, x1]
      else return(NULL)
    })
    combine_info_DA <- do.call(base::cbind, tmp1)
    tmp1 <- lapply(col_ord, function(x) {
      combine_info_DE[, grep(sprintf("^%s\\.", x), colnames(combine_info_DE))]
    })
    combine_info_DE <- do.call(base::cbind, tmp1)
  }
  ms_tab <- base::cbind(label_info, combine_info_DA, combine_info_DE, 
                        add_info)
  rownames(ms_tab) <- ms_tab$originalID_label
  return(ms_tab)
}


############################
#    COMMAND LINE OPTIONS
############################
option_list <- list(
  make_option(c("-d", "--project_main_dir"), type="character", default="./",
              help="Main directory for the project [default = %default]"),
  
  make_option(c("-p", "--project_name"), type="character", default="Driver_Inference",
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
#source("bin/Adjusted_functions/NetBID.analysis.dir.create_nc.R")

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

#source("bin/Adjusted_functions/merge_TF_SIG_NC.network.R")

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

# ---- REMOVE BIOMART ----
analysis.par$transfer_tab <- NULL

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

out_activity <- file.path(
  analysis.par$out.dir.DATA,
  paste0(analysis.par$project.name, "_activity_matrix.csv")
)

out2excel(analysis.par$final_ms_tab, out.xlsx = out_file)

write.csv(
  as.data.frame(exprs(analysis.par$merge.ac.eset)),
  file = out_activity,
  row.names = TRUE,
  quote = FALSE
)

NetBID.saveRData(analysis.par = analysis.par, step = "ms-tab")

cat("Analysis completed successfully.\n")
cat("Output saved to:", out_file, "\n")