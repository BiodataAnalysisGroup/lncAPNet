process DRIVER_INFERENCES {

    container "vasileioubill/netbid2:latest"

    input:
    tuple val(meta), path(netbid_dir), path(sjaracne_dir)

    output:
    tuple val(meta), path("Driver_Inference"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy the entire netbid_dir structure to Driver_Inference
    mkdir Driver_Inference
    cp -r ${netbid_dir}/* Driver_Inference/
    mv Driver_Inference/SJAR/NetBID2_Project Driver_Inference/SJAR/Driver_Inference
    
    # Copy sjaracne results into the existing SJAR/NetBID2_Project directory
    cp -r ${sjaracne_dir}/* Driver_Inference/SJAR/Driver_Inference/
    
    Driver_inference.R \\
        --project_main_dir ./ \\
        --project_name Driver_Inference
    """
}