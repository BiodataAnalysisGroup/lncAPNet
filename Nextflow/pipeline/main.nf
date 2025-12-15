#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NETWORK_RECONSTRUCTION } from './modules/local/network_reconstruction/main'

workflow {
    ch_eset = channel.of([
        [ id: 'eset' ],
        file(params.input, checkIfExists: true)
    ])

    NETWORK_RECONSTRUCTION( ch_eset, params.iqr )


}
