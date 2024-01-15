#!/usr/bin/env nextflow

process CUTADAPT_Trim1_or {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fastq.gz'), emit: trim_out
    //tuple val(sample), file('*_I7Trim.fastq.gz')
    script:
    """
    cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -A GTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 46 --pair-filter=first --overlap 18 -o ${sample}_i7Trim_R1.fastq.gz -p ${sample}_i7Trim_R2.fastq.gz ${pe_reads[0]} ${pe_reads[1]} 
    """
}

process CUTADAPT_Trim4 {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
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
    publishDir "${params.outdir}/${index_step}_cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
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
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    val (sample)
    tuple path(r1), path(r2)

    output:
    tuple val(sample), file('*.fastq.gz'), emit: trim_out
    //tuple val(sample), file('*_I7Trim.fastq.gz')
    script:
    """
    cutadapt -a T{16} --minimum-length 46 --pair-filter=first --overlap 15 -o ${sample}_R1_Ttrim.fastq.gz -p ${sample}_Ttrim_R2.fastq.gz ${r1} ${r2} 
    """
}


process CUTADAPT_QSP {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trimmed.fastq.gz')
    script:
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[0]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample}_trimmed.fastq.gz -
    """
}

process CUTADAPT_QSP_2 {

    container './containers/cutadapt:4.5--py39hf95cd2a_0.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trimmed.fastq.gz')
    script:
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[1]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample}_trimmed.fastq.gz -
    """
}
