#!/usr/bin/env nextflow

process plot_reads {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true 
    //label 'high'
    module 'R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1'

    input:
    path (files)
    output:
    file ("*")

    script:
    """
    CountReads.r ${files} 
    """
    }
