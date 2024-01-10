#!/usr/bin/env nextflow

process Decode {
    publishDir "${params.outdir}/${index_step}_demux", mode: 'copy', overwrite: true 
    //label 'high'
    module 'PyTorch/1.12.0-foss-2022a-CUDA-11.7.0'
    tag "${sample}" 
    
    input:
    val (index_step)
    path (B)
    val (N)
    val (Nthr)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file("*.fastq.gz")
    
    script:
    """
    Demux.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${N} -M ${params.M} -B ${B} -P ${params.P} --Ntriage ${params.Ntriage} --Nthresh ${Nthr} --out-prefix ${sample} 
    """
}



