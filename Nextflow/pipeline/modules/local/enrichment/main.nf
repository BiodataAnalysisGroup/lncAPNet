process ENRICHMENT {

    container "vasileioubill/pasnet:latest"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path(gmt_files)  // This receives the collected GMT files as a list
    tuple val(meta), path(driver_dir)

    output:
    tuple val(meta), path("Enrichment"), emit: results
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Sort GMT files to ensure consistent ordering
    def gmt_list = gmt_files instanceof List ? gmt_files : [gmt_files]
    def gmt_GO = gmt_list.find { file -> file.name.contains('GO_Biological_Process') }
    def gmt_Reactome = gmt_list.find { file -> file.name.contains('Reactome') }
    def gmt_KEGG = gmt_list.find { file -> file.name.contains('KEGG') }
    def gmt_WikiPathway = gmt_list.find { file -> file.name.contains('WikiPathway') }
    """
    mkdir -p Enrichment/PASNet/Input/GO
    mkdir -p Enrichment/PASNet/Input/REST

    # Copy driver inference results
    cp ${driver_dir}/Driver_output/${driver_dir}/DATA/Driver_Inference_ms_tab.xlsx Enrichment/
    cp ${driver_dir}/Driver_output/${driver_dir}/DATA/Driver_Inference_metadata.csv Enrichment/
    cp ${driver_dir}/Driver_output/${driver_dir}/DATA/Driver_Inference_activity_matrix.csv Enrichment/

    # Run pathway enrichment
    Pathway_Enrichment.py \\
        --ms_tab Enrichment/Driver_Inference_ms_tab.xlsx \\
        --metadata Enrichment/Driver_Inference_metadata.csv \\
        --gmt_GO ${gmt_GO} \\
        --gmt_Reactome ${gmt_Reactome} \\
        --gmt_KEGG ${gmt_KEGG} \\
        --gmt_Wikipathway ${gmt_WikiPathway} \\
        --outdir Enrichment
    """
}