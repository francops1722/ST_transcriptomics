
process FeatureCounts {
    container './containers/subread:2.0.1--hed695b0_0.sif'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    
    input:
    val(index_step)
    tuple val(sample), path(bam_file)
    path gtf_file
    path idx

    output:
    path("*.count"), emit: counts
    path("*.summary"), emit: logs

    script:
    """
    featureCounts -s 1 -T 10 -t exon -g gene_id -a ${gtf_file}  -o ${sample}_uniquev5.count ${bam_file} 
    """
}

process FeatureCounts_BAM {
    container './containers/subread:2.0.1--hed695b0_0.sif'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    path gtf_file
    path idx

    output:
    tuple val(sample), path("*.bam"), emit: FC_bam
    path("*.summary"), emit: log_files

    script:
    """
    featureCounts -s 1 -T 4 -t exon -g gene_id -a ${gtf_file}  -o ${sample} -R BAM ${bam_file} 
    """
}

process index_bam_FC {

    container './containers/samtools:1.16.sif'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true
    tag "${sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    output:
    path("*.bai")

    script:
    """
    samtools index ${bam_file}
    """
}

process sort_bam_FC {

    container './containers/samtools:1.16.sif'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true
    tag "${sample}"
    
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

process merge_featureCounts {
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy'

    input:
    val (index_step)
    file input_files 

    output:
    file 'merged_gene_counts.txt'

    script:
    """
    FC_merge.py -o merged_gene_counts.txt -i $input_files
    """
}

process merge_Counts {
    publishDir "${params.outdir}/${index_step}_Counts_summary", mode: 'copy'
    //module 'pandas/1.1.2-foss-2020a-Python-3.8.2' only use when working in cluster joltik 
    module 'PyTorch/1.12.0-foss-2022a-CUDA-11.7.0'

    input:
    val (index_step)
    file input_files 

    output:
    file 'gene_counts.csv'
    file 'summary_barcodes.csv'
    val 'ready', emit: mock

    script:
    """
    Counts_merge.py -i $input_files -Bar summary_barcodes.csv -Count gene_counts.csv
    """
}

