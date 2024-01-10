#!/usr/bin/env nextflow

process Subsample_Seqtk {
    container './containers/seqtk:1.3--hed695b0_2.sif'
    publishDir "${params.outdir}/${index_step}_Subsampled", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
    val(sample)
    tuple path(r1), path(r2)

    output:
    tuple path("*R1_subs.fastq.gz"), path("*R2_subs.fastq.gz")

    script:
    """
    seqtk sample -s 100  ${r1} ${params.Subs} | gzip > ${sample}_R1_subs.fastq.gz;
	seqtk sample -s 100  ${r2} ${params.Subs} | gzip > ${sample}_R2_subs.fastq.gz;
    """
}

process Subsample_Seqtk_or {
    container './containers/seqtk:1.3--hed695b0_2.sif'
    publishDir "${params.outdir}/${index_step}_Subsampled", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_subs.fastq.gz')

    script:
    """
    seqtk sample -s 100  ${pe_reads[0]} ${params.Subs} | gzip > ${sample}_R1_subs.fastq.gz;
	seqtk sample -s 100  ${pe_reads[1]} ${params.Subs} | gzip > ${sample}_R2_subs.fastq.gz;
    """
}


process Subsample_Seqtk_SE {
    container './containers/seqtk:1.3--hed695b0_2.sif'
    publishDir "${params.outdir}/Subsampled", mode: 'copy', overwrite: true

    input:
    path(se_reads)

    output:
    file('*_subs.fastq.gz')

    script:
    """
    seqtk sample -s 100  ${se_reads} ${params.Subs} | gzip > _subs.fastq.gz
    """
}

