process PASNET {

    container "vasileioubill/pasnet:latest"

    input:
    tuple val(meta), path(enrichment)

    output:
    tuple val(meta), path("PASNet"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy the entire netbid_dir structure to Driver_Inference
    mkdir PASNet
    cp -r ${enrichment}/PASNet/* PASNet/
    mkdir -p PASNet/Output
    
    PASNet_run.py \\
      --train-data PASNet/Input/GO/Training.xlsx \
      --val-data PASNet/Input/GO/Validation.xlsx \
      --pathway-mask PASNet/Input/GO/pt_fixed.xlsx \
      --in-nodes 1756 \
      --pathway-nodes 140 \
      --hidden-nodes 140 \
      --output-dir results \
      --plot-roc \
      --shap-analysis
    """
}