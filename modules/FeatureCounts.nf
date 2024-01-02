
process FeatureCounts_BAM {
    container './containers/subread:2.0.1--hed695b0_0.sif'
    publishDir "${params.outdir}/Counts", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    
    input:
    tuple val(sample), path(bam_file)
    path gtf_file
    path idx

    output:
    tuple val(sample), path("*.bam"), emit: FC_bam
    path("*.summary"), emit: logs

    script:
    """
    featureCounts -s 1 -T 4 -t exon -g gene_id -a ${gtf_file}  -o ${sample} -R BAM ${bam_file} 
    """
}

process index_bam_FC {

    container './containers/samtools:1.16.sif'
    publishDir "${params.outdir}/Counts", mode: 'copy', overwrite: true
    tag "${sample}"
    
    input:
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
    publishDir "${params.outdir}/Counts", mode: 'copy', overwrite: true
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam_file)
    
    output:
    tuple val(sample), path("*.bam"), emit: sorted_bam

    script:
    """
    samtools sort ${bam_file} -o ${sample}_assigned_sorted.bam
    """
}

process merge_featureCounts {
    publishDir "${params.outdir}/Counts", mode: 'copy'

    input:
    file input_files 

    output:
    file 'merged_gene_counts.txt'

    script:
    """
    FC_merge.py -o merged_gene_counts.txt -i $input_files
    """
}



