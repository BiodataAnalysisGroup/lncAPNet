process ENRICHMENT {

    container "vasileioubill/pasnet:latest"

    input:
    tuple val(meta), path(driver_dir), 

    output:
    tuple val(meta), path("Enrichment"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p PASNet/Input/GO
    mkdir -p PASNet/Input/REST
    mkdir -p Enrichment/
    mkdir -p NetBID2_output/

    cp ${driver_dir}/Driver_output/Data/*.xlsx NetBID2_output/

    Pathway_Enrichment.py \\
        --project_dir ./ \\
        --project_name Enrichment
    """
}