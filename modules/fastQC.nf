#!/usr/bin/env nextflow

process FASTQC { 
    // DIRECTIVES: set the docker container, the directory to output to, and a tag to follow along which sample is currently being processed
    container './containers/fastqc:0.11.9--0.sif'
    publishDir "${params.outdir}/fastqc/${step}", mode: 'copy', overwrite: true
    tag "${step}"

    input:
    val (step)
    tuple path(r1), path(r2)

    output:
    path('*_fastqc.{zip,html}')

    script:
    """
    fastqc ${r1} ${r2}
    """
}

process MULTIQC {
    // DIRECTIVES: set the docker container, the directory to output to, and a tag
    container './containers/multiqc:latest.sif'
    publishDir "${params.outdir}/multiqc/${step}/", mode: 'copy', overwrite: true
    tag "${step}"

    input:
    val(step)
    path(inputfiles)

    output:
    file('*_multiqc_report.html')

    script:
    """
    multiqc --filename ${step}_multiqc_report .
    """
}


// combine the two processes into a subworkflow
workflow check_QC {
    take:
        step
        reads_ch

    main:
        FASTQC(step, reads_ch)
        
        input_multiqc = FASTQC.out.collect()
        MULTIQC(step, input_multiqc)
}
