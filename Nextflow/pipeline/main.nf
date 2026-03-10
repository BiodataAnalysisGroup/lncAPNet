#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NETWORK_RECONSTRUCTION } from './modules/local/network_reconstruction/main'
include { SJARACNE               } from './modules/local/sjaracne/main'
include { DRIVER_INFERENCES      } from './modules/local/driver_inferences/main'
include { ENRICHMENT             } from './modules/local/enrichment/main'

workflow {
    ch_eset = channel.of([
        [ id: 'eset' ],
        file(params.input, checkIfExists: true)
    ])

    ch_gene_info = channel.of(file(params.gene_info, checkIfExists: true))

    ch_gmt_files = channel.fromPath('bin/Enrichment/*.gmt', checkIfExists: true).collect()

    NETWORK_RECONSTRUCTION( ch_eset, ch_gene_info.first(), params.iqr )

    SJARACNE( NETWORK_RECONSTRUCTION.out.results )
    
    DRIVER_INFERENCES( NETWORK_RECONSTRUCTION.out.results.join(SJARACNE.out.results) )

    ENRICHMENT( DRIVER_INFERENCES.out.results )

}