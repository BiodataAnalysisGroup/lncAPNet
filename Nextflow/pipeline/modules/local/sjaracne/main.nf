process SJARACNE {

    container "vasileioubill/sjaracne:latest"

    input:
    tuple val(meta), path(netbid_dir)

    output:
    tuple val(meta), path("*"), emit: results

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    SJARACNe_run.sh \\
        --project_dir ${netbid_dir} \\
        --project_name ${prefix}
    """
}
