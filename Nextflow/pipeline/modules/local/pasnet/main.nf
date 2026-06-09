process PASNET {

    container "vasileioubill/pasnet:latest"

    input:
    val(comparison)
    tuple val(meta), path(enrichment)

    output:
    tuple val(meta), path("PASNet"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir PASNet
    cp -r ${enrichment}/PASNet/* PASNet/
    
    PASNet_run.py \
      --train-data PASNet/Input/GO/Training.xlsx \
      --val-data PASNet/Input/GO/Validation.xlsx \
      --val-data-grid PASNet/Input/GO/Validation_grinding.xlsx \
      --pathway-mask PASNet/Input/GO/pt_fixed.xlsx \
      --comparison ${comparison} \
      --outdir PASNet
    """
}
