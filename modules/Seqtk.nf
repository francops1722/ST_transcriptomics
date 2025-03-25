#!/usr/bin/env nextflow

process Subsample_Seqtk {
    container './containers/seqtk.sif'
    label 'med'
    publishDir "${params.outdir}/${index_step}_Subsampled", mode: 'copy', overwrite: true
    tag "${sample.sample}"

    input:
    val (index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fastq.gz')

    script:
    """
    seqtk sample -s 100  ${pe_reads[0]} ${params.Subs} | gzip > ${sample.sample}_subs_R1.fastq.gz;
	seqtk sample -s 100  ${pe_reads[1]} ${params.Subs} | gzip > ${sample.sample}_subs_R2.fastq.gz;
    """
}


process Subsample_Seqtk_or {
    container './containers/seqtk.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_Subsampled", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fastq.gz')

    script:
    """
    seqtk sample -s 100  ${pe_reads[0]} ${params.Subs} | gzip > ${sample}_subs_R1.fastq.gz;
	seqtk sample -s 100  ${pe_reads[1]} ${params.Subs} | gzip > ${sample}_subs_R2.fastq.gz;
    """
}


process Subsample_Seqtk_SE {
    container './containers/seqtk.sif'
    //label 'med'
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

