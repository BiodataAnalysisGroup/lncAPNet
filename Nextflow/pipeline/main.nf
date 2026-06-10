#!/usr/bin/env python
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NETWORK_RECONSTRUCTION } from './modules/local/network_reconstruction/main'
include { SJARACNE               } from './modules/local/sjaracne/main'
include { DRIVER_INFERENCES      } from './modules/local/driver_inferences/main'
include { ENRICHMENT             } from './modules/local/enrichment/main'
include { PASNET                 } from './modules/local/pasnet/main'

workflow {
    ch_eset = channel.of([
        [ id: 'eset' ],
        file(params.input, checkIfExists: true)
    ])

    ch_gene_info = channel.of(file(params.gene_info, checkIfExists: true))

    NETWORK_RECONSTRUCTION( ch_eset, ch_gene_info.first(), params.iqr )

    SJARACNE( NETWORK_RECONSTRUCTION.out.results )
    
    // Group Patameters values
    ch_group0 = channel.of(params.group0)
    ch_group1 = channel.of(params.group1)
    ch_comparison = channel.of(params.comparison)

    DRIVER_INFERENCES(ch_group0, ch_group1, ch_comparison, NETWORK_RECONSTRUCTION.out.results.join(SJARACNE.out.results) )

    // gmt pathway-genes files for pathway Enrichment
    ch_gmt_files = channel.fromPath('bin/Enrichment/*.gmt', checkIfExists: true).collect()

    ENRICHMENT(ch_gmt_files, ch_group0, ch_group1, ch_comparison, DRIVER_INFERENCES.out.results)

    // Trained Model saved as pickle file

    PASNET(ch_comparison, ENRICHMENT.out.results)

}