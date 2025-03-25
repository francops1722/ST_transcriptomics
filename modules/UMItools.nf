#!/usr/bin/env nextflow
process UMI_extract {
    container './containers/umi_tools.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_UMI_extracted", mode: 'copy', overwrite: true 
    tag "${sample.sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file('*_{R1,R2}_UMIextracted.fastq.gz'), emit: fastq
    path('logs/*.log'), emit: log_files
    

    script:
    if ("$params.UMItype" == 'semi')
        """
        mkdir logs
        umi_tools extract --stdin ${pe_reads[0]} --stdout ${sample.sample}_R1_UMIextracted.fastq.gz -L logs/${sample.sample}_extraction.log --extract-method=string --bc-pattern NNNNNNNNN --read2-in ${pe_reads[1]} --read2-out ${sample.sample}_R2_UMIextracted.fastq.gz 
        """
    else if ("$params.UMItype" == 'random')
        """
        mkdir logs
        umi_tools extract --stdin ${pe_reads[0]} --stdout ${sample.sample}_R1_UMIextracted.fastq.gz -L logs/${sample.sample}_extraction.log --extract-method=string --bc-pattern NNNNNNNNNN --read2-in ${pe_reads[1]} --read2-out ${sample.sample}_R2_UMIextracted.fastq.gz 
        """
    else
        error "Invalid UMI pattern (use only: semi or random)"
}

//to do like QSP - where read 2 has the UMi and barcode and read 1 is for mapping
process UMI_QSP {
    container './containers/umi_tools.sif'
    publishDir "${params.outdir}/${index_step}_UMI_extracted", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file('*_UMIextracted.fastq.gz'), emit: fastq
    path('logs/*_log'), emit: logs
    

    script:
    """
    mkdir logs
    umi_tools extract --stdin ${pe_reads[0]} --stdout ${sample}_R1_UMIextracted.fastq.gz -L logs/${sample}_extraction_log --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN --read2-in ${pe_reads[1]} --read2-out ${sample}_R2_UMIextracted.fastq.gz 
    """   
}

//Auxiliary function similar to QSP but read 1 has the UMi and barcode and read 2 is for mapping
process UMI_QSP_2 {
    container './containers/umi_tools.sif'
    publishDir "${params.outdir}/${index_step}_UMI_extracted", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file('*_UMIextracted.fastq.gz'), emit: fastq
    path('logs/*_log'), emit: logs
    

    script:
    """
    mkdir logs
    umi_tools extract --stdin ${pe_reads[1]} --stdout ${sample}_R2_UMIextracted.fastq.gz -L logs/${sample}_extraction_log --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN --read2-in ${pe_reads[0]} --read2-out ${sample}_R1_UMIextracted.fastq.gz 
    """   
}
process UMI_count {

    container './containers/umi_tools.sif'
    publishDir "${params.outdir}/${index_step}_UMI_Counts", mode: 'copy', overwrite: true 
    tag "${sample.sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    tuple val(sample), path(index)
    
    
    output:
    path('*.tsv.gz'), emit: umi_counts
    path('logs_umi/*.{txt,log}'), emit: log_files

    script:
    """
    mkdir logs_umi
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bam_file} -S ${sample.sample}_counts.tsv.gz --log=logs_umi/${sample.sample}_umiCounts.log
    """

}

// process UMI_dedup_basic {

//     container './containers/umi_tools.sif'
//     publishDir "${params.outdir}/${index_step}_UMI_dedup", mode: 'copy', overwrite: true 
//     tag "${sample}"
    
//     input:
//     val (index_step)
//     tuple val(sample), path(bam_file)
//     //tuple val(sample), path(index)
    
//     output:
//     //path("*")
//     tuple val(sample), path('*.bam'), emit: dedup_files
//     path('logs_dedup/*.{txt,log}'), emit: logs

//     script:
//     """
//     mkdir logs_dedup
//     umi_tools dedup -I ${bam_file} -S ${sample}_dedup.bam --multimapping-detection-method=NH --output-stats=logs_dedup/${sample}_deduplicated.txt --log=logs_dedup/${sample}_deduplication.log
//     """

// }

process UMI_dedup_basic {

    container './containers/umi_tools.sif'
    publishDir "${params.outdir}/${index_step}_Star", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file), path(index)
    //tuple val(sample), path(index)

    output:
    //path("*")
    tuple val(sample), path('*.bam'), emit: dedup_files
    path('logs_dedup/*.{txt,log}'), emit: logs

    script:
    """
    mkdir logs_dedup
    umi_tools dedup -I ${bam_file} -S ${sample}_dedup.bam --multimapping-detection-method=NH --output-stats=logs_dedup/${sample}_deduplicated.txt --log=logs_dedup/${sample}_deduplication.log
    """

}