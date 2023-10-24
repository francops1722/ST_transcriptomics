#!/usr/bin/env nextflow

process Subsample_Seqtk {
    container './containers/seqtk:1.4--he4a0461_1.sif'
    publishDir "${params.outdir}/Subsampled", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
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

