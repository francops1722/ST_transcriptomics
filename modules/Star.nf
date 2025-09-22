#!/usr/bin/env nextflow

process Star_Align_R2 {
    
    container './containers/star.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true 
    label "high"
    tag "${sample}"
    
    
    input:
    val(index_step)
    tuple val(sample), path(se_reads)
    path genome
    
    output:
    path("*_Log.final.out"), emit: log_files
    tuple val(sample), path("*.bam"), emit: align_bam

    script:
    """
    STAR --readFilesCommand zcat --genomeDir ${genome} --readFilesIn ${se_reads} --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.5 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitOutSJcollapsed 5000000 --limitIObufferSize 200000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}_ --limitBAMsortRAM 2000000000
    """
}

process Star_Align_R2v2 {
    
    container './containers/star.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true 
    tag "${sample.sample}"
    label "high"
    
    input:
    val (index_step)
    tuple val(sample), path(se_reads)
    path genome
    
    output:
    path("*_Log.final.out"), emit: log_files
    tuple val(sample), path("*.bam"), emit: align_bam

    script:
    """
    STAR --readFilesCommand zcat --genomeDir ${genome} --readFilesIn ${se_reads} --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitOutSJcollapsed 5000000 --limitIObufferSize 200000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample.sample}_ --limitBAMsortRAM 2000000000
    """
}

process Star_Align_QSP {
    
    container './containers/star.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true 
    tag "${sample}"
    label "high"
    
    input:
    val(index_step)
    tuple val(sample), path(se_reads)
    path genome
    
    output:
    path("*_Log.final.out"), emit: log_files
    tuple val(sample), path("*sortedByCoord.out.bam"), emit: align_bam
    
    script:
    """
    STAR --readFilesCommand zcat --genomeDir ${genome} --readFilesIn ${se_reads} --outFilterType BySJout --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitOutSJcollapsed 5000000 --limitIObufferSize 200000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}_ --limitBAMsortRAM 2000000000 
    """
}


process Star_Align_QSP_test {
    
    container './containers/star.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true 
    tag "${sample}"
    label "high"
    
    input:
    val(index_step)
    tuple val(sample), path(se_reads)
    path genome
    
    output:
    path("*_Log.final.out"), emit: log_files
    tuple val(sample), path("*.bam"), emit: align_bam
    
    script:
    """
    STAR \
    --twopassMode Basic \
    --readFilesCommand zcat \
    --genomeDir ${genome} \
    --readFilesIn ${se_reads} \
    --outFilterType BySJout \
    --outFilterMultimapNmax 200 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.6 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --limitOutSJcollapsed 5000000 \
    --limitIObufferSize 200000000 \
    --outSAMattributes NH HI NM MD \
    --outSAMunmapped Within \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${sample}_ \
    --limitBAMsortRAM 2000000000

"""
}


process index_bam {
    

    container './containers/samtools.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true
    tag "${sample.sample}"
    label "low"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    output:
    tuple val(sample), path("*.bai")

    script:
    """
    samtools index ${bam_file}
    """
}

process index_bam_QSP {
    container './containers/samtools.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true
    tag "${sample}"
    label "low"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    output:
    tuple val(sample), path("*.bai")
    //val ("ready"), emit: mock 
    
    script:
    """
    samtools index ${bam_file}
    """
}

process sort_bam {

    container './containers/samtools.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true
    tag "${sample.sample}"
    label "low"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    
    output:
    tuple val(sample), path("*.bam"), emit: sorted_bam

    script:
    """
    samtools sort ${bam_file} -o ${sample}_assigned_sorted.bam
    """
}

process sort_bam_QSP {

    container './containers/samtools.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true
    tag "${sample}"
    label "low"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    
    output:
    tuple val(sample), path("*.bam"), emit: sorted_bam

    script:
    """
    samtools sort ${bam_file} -o ${sample}_assigned_sorted.bam
    """
}




