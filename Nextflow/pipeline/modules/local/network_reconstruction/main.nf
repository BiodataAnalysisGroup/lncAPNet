process NETWORK_RECONSTRUCTION {

    container "vasileioubill/netbid2:latest"

    input:
    tuple val(meta), path(rds)
    val(iqr)

    output:
    tuple val(meta), path("NetBID2_Project"), emit: results

    script:
    """
    Network_reconstruction.R \\
        --eset ${rds} \\
        --project_name NetBID2_Project \\
        --iqr ${iqr}
    """
}
