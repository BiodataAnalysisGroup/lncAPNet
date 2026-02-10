process NETWORK_RECONSTRUCTION {

    container "vasileioubill/netbid2:latest"

    input:
    tuple val(meta), path(rds)
    path(gene_info)
    val(iqr)

    output:
    tuple val(meta), path("NetBID2_Project"), emit: results

    script:
    """
    Network_reconstruction.R \\
        --eset ${rds} \\
        --gene_info ${gene_info} \\
        --project_name NetBID2_Project \\
        --iqr ${iqr}
    """
}
