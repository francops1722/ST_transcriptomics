
process FeatureCounts {
    container './containers/subread.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true   
    tag "${sample}" 
    
    input:
    val(index_step)
    tuple val(sample), path(bam_file)
    path gtf_file
    //tuple val(sample), path(bam_file_or), path(index)
    //tuple val(sample), path(index)

    output:
    path("*.count"), emit: counts
    path("*.summary"), emit: logs

    script:
    """
    featureCounts -s 1 -T 10 -t exon -g gene_id -a ${gtf_file}  -o ${sample}_uniquev5.count ${bam_file} 
    """
}

process FeatureCounts_transcript_id {
    container './containers/subread.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true   
    tag "${sample}" 
    
    input:
    val(index_step)
    tuple val(sample), path(bam_file)
    path gtf_file
    //tuple val(sample), path(bam_file_or), path(index)
    //tuple val(sample), path(index)

    output:
    path("*.count"), emit: counts
    path("*.summary"), emit: logs

    script:
    """
    featureCounts -s 1 -T 10 -t exon -g transcript_id -a ${gtf_file}  -o ${sample}_uniquev5.count ${bam_file} 
    """
}


process FeatureCounts_BAM {
    container './containers/subread.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true 
    tag "${sample.sample}"
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    path gtf_file
    tuple val(sample), path(idx)

    output:
    tuple val(sample), path("*.bam"), emit: FC_bam
    path("*.summary"), emit: log_files

    script:
    """
    featureCounts -s 1 -T 4 -t exon -g gene_id -a ${gtf_file}  -o ${sample.sample} -R BAM ${bam_file} 
    """
}

process index_bam_FC {

    container './containers/samtools.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true
    tag "${sample.sample}"
    
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

process sort_bam_FC {
    container './containers/samtools.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_FCounts", mode: 'copy', overwrite: true
    tag "${sample.sample}" 
    
    input:
    val (index_step)
    tuple val(sample), path(bam_file)
    
    output:
    tuple val(sample), path("*.bam"), emit: sorted_bam

    script:
    """
    samtools sort ${bam_file} -o ${sample.sample}_assigned_sorted.bam
    """
}

process merge_featureCounts {
    publishDir "${params.outdir}/${index_step}_Counts_summary", mode: 'copy'
    label 'low'

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

process merge_TPMs {
    publishDir "${params.outdir}/${index_step}_Counts_summary", mode: 'copy'
    label 'low'

    input:
    val (index_step)
    file input_files

    output:
    file 'merged_salmon_TPM.txt'


    script:
    """
    merge_salmon_quant.py -i $input_files
    """
}

process merge_Counts {
    publishDir "${params.outdir}/${index_step}_Counts_summary", mode: 'copy'
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



process table_reads {
    publishDir "${params.outdir}/${index_step}_Counts_summary", mode: 'copy', overwrite: true 
    module 'R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1'

    input:
    val (index_step)
    path files
    val ready
    output:
    file('*')
    
    script:
    """
    MakeReadsTable.r ${files}
    """
    }
