#!/usr/bin/env nextflow

process CUTADAPT_Trim1_or {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_i7Trim.fastq.gz'), emit: trim_out
    //tuple val(sample), file('*_I7Trim.fastq.gz')
    script:
    """
    cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -A GTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 46 --pair-filter=first --overlap 18 -o ${sample}_R1_i7Trim.fastq.gz -p ${sample}_R2_i7Trim.fastq.gz ${pe_reads[0]} ${pe_reads[1]} 
    """
}

process CUTADAPT_Trim4 {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 15 -q 20 -O 20 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA ${pe_reads[1]}| cutadapt -m 15 "-a polyA=A{18}" -o ${sample}_Trim4_R2.fastq.gz -;
    """
}

process CUTADAPT_Trim4_1 {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 10 -q 20 "-a polyA=A{18}" -o ${sample}_Trim4_R2.fastq.gz ${pe_reads[1]}
    """
}

process CUTADAPT_Trim5 {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 15 -q 20 -O 20 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA ${pe_reads[1]}| cutadapt -m 15 "-a polyA=A{18}" -o ${sample}_Trim4_R2.fastq.gz -;
    """
}

process CUTADAPT_oligodT {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (sample)
    tuple path(r1), path(r2)

    output:
    tuple val(sample), file('*_Ttrim.fastq.gz'), emit: trim_out
    //tuple val(sample), file('*_I7Trim.fastq.gz')
    script:
    """
    cutadapt -a T{16} --minimum-length 46 --pair-filter=first --overlap 15 -o ${sample}_R1_Ttrim.fastq.gz -p ${sample}_R2_Ttrim.fastq.gz ${r1} ${r2} 
    """
}

