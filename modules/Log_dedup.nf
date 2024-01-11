#!/usr/bin/env nextflow

process Dedup_log {
    publishDir "${params.outdir}/${index_step}_UMI_Counts", mode: 'copy', overwrite: true 
    tag "${sample}" 
    
    input:
    val(index_step)
    path(log_files)
    
    script:
    """
    for sample in ${log_files}; 
    do inputreads=`grep -i "input reads" \$sample | awk 'NF>1{print \$NF}'`; 
    outputreads=`grep -i "reads counted" \$sample | awk 'NF>1{print \$NF}'`;
    echo \$sample \$inputreads \$outputreads | cat >> ${params.outdir}/${index_step}_UMI_Counts/logs_umi/UMI_dedup.txt; done
    """
}



