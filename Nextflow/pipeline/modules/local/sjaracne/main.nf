process SJARACNE {

    container "vasileioubill/sjaracne:latest"

    input:
    tuple val(meta), path(netbid_dir)

    output:
    tuple val(meta), path("SJARACNe"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p SJARACNe
    cp ${netbid_dir}/SJAR/${netbid_dir}/* SJARACNe/

    SJARACNe_run.sh \\
        --project_dir ./ \\
        --project_name SJARACNe
    """
}