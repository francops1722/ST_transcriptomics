#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp(""" Usage: nextflow run main.nf 
   
   Pipeline to analyze RNAseq data obtained from
   our custom Spatial transcriptomics platform
   """)
   System.exit(0)
}
// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

//main workflow
include { ST_workflow } from './workflows/ST_main'
include { QSP_workflow } from './workflows/QSP2_main'
include { QSPLike_workflow } from './workflows/QSPLike_main'

workflow{
    if (params.flavour=="ST"){
        ST_workflow ()
    } else if (params.flavour=="QSP") {
        QSP_workflow()
    } else {
        QSPLike_workflow()
    }
    
}

