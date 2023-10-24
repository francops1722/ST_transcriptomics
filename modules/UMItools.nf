#!/usr/bin/env nextflow
process UMI_extract {
    container './containers/umitools:latest.sif'
    publishDir "${params.outdir}/UMI_extracted", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file('*_UMIextracted.fastq.gz'), emit: fastq
    path('logs/*_log'), emit: logs
    

    script:
    if ("$params.UMItype" == 'semi')
        """
        mkdir logs
        umi_tools extract --stdin ${pe_reads[0]} --stdout ${sample}_R1_UMIextracted.fastq.gz -L logs/${sample}_extraction_log --extract-method=string --bc-pattern NNNNNNNNN --read2-in ${pe_reads[1]} --read2-out ${sample}_R2_UMIextracted.fastq.gz 
        """
    else if ("$params.UMItype" == 'random')
        """
        mkdir logs
        umi_tools extract --stdin ${pe_reads[0]} --stdout ${sample}_R1_UMIextracted.fastq.gz -L logs/${sample}_extraction_log --extract-method=string --bc-pattern NNNNNNNNNN --read2-in ${pe_reads[1]} --read2-out ${sample}_R2_UMIextracted.fastq.gz 
        """
    else
        error "Invalid UMI pattern (use only: semi or random)"
}

process UMI_count {

    container './containers/umitools:latest.sif'
    publishDir "${params.outdir}/Counts", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam_file)
    path idx
    
    output:
    path('*')
    
    script:
    """
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bam_file} -S ${sample}_counts.tsv.gz
    """

}

process UMI_dedup_basic {

    container './containers/umitools:latest.sif'
    publishDir "${params.outdir}/UMI_dedup", mode: 'copy', overwrite: true 
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam_file)
    path index_file
    
    output:
    //path("*")
    path('*.bam'), emit: dedup_files
    path('logs_dedup/*.{txt,log}'), emit: logs

    script:
    """
    mkdir logs_dedup
    umi_tools dedup -I ${bam_file} -S ${sample}_dedup.bam --multimapping-detection-method=NH --output-stats=logs_dedup/${sample}_deduplicated.txt --log=logs_dedup/${sample}_deduplication.log
    """

}
