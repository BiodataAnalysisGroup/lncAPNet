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
