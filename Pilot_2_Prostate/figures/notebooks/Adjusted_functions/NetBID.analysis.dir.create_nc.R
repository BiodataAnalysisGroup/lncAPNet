
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
