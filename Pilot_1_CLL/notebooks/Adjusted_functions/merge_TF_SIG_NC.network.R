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
