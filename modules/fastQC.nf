#!/usr/bin/env nextflow

process FASTQC { 
    // DIRECTIVES: set the docker container, the directory to output to, and a tag to follow along which sample is currently being processed
    container './containers/fastqc.sif'
    label 'low'
    publishDir "${params.outdir}/fastqc/${step}", mode: 'copy', overwrite: true
    tag "${step}"

    input:
    val(step)
    tuple val(sample), path(reads)
    
    output:
    path('*_fastqc.{zip,html}')

    script:
    """
    fastqc ${reads}
    """
}


process MULTIQC {
    // DIRECTIVES: set the docker container, the directory to output to, and a tag
    container './containers/multiqc.sif'
    label 'low'
    publishDir "${params.outdir}/multiqc/${step}/", mode: 'copy', overwrite: true
    tag "${step}"

    input:
    val(step)
    path(inputfiles)

    output:
    file('*')

    script:
    """
    multiqc --filename ${step}_multiqc_report .
    """
}


