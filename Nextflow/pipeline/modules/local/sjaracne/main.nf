process SJARACNE {

    container "vasileioubill/sjaracne:latest"

    input:
    tuple val(meta), path(netbid_dir)

    output:
    tuple val(meta), path("*"), emit: results

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${netbid_dir}/eset/SJAR/eset/
    cp ${netbid_dir}/SJAR/${netbid_dir}/* ${netbid_dir}/eset/SJAR/eset/
    
    SJARACNe_run.sh \\
        --project_dir ${netbid_dir} \\
        --project_name ${prefix}
    """
}
