#!/usr/bin/env nextflow

process main_trim {
    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_4", mode: 'copy', overwrite: true

    input:
    val(index_step)
    tuple val (sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz'), emit: trim_out

    script:
    """
    cutadapt -m 15 -q 20 -O 20 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
    -A polyA=A{18} \
    -o ${sample.sample}_Trim_R1.fastq.gz -p ${sample.sample}_Trim_R2.fastq.gz ${pe_reads[0]} ${pe_reads[1]}
    """
}


process CUTADAPT_Trim1_or {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fastq.gz'), emit: trim_out
    script:
    """
    cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG -A GTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 46 --pair-filter=first --overlap 18 -o ${sample}_i7Trim_R1.fastq.gz -p ${sample}_i7Trim_R2.fastq.gz ${pe_reads[0]} ${pe_reads[1]} 
    """
}

process CUTADAPT_Trim2 {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_2", mode: 'copy', overwrite: true
    tag "${sample.sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[1]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample.sample}_R2_trim2.fastq.gz -
    """
}

process CUTADAPT_Trim4_1 {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 10 -q 20 "-a polyA=A{18}" -o ${sample}_Trim4_R2.fastq.gz ${pe_reads[1]}
    """
}

process CUTADAPT_Trim5 {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/cutadapt_4", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)
    

    output:
    tuple val(sample), file('*.fastq.gz')
    
    script:
    """
     cutadapt -m 15 -q 20 -O 20 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA ${pe_reads[1]}| cutadapt -m 15 "-a polyA=A{18}" -o ${sample}_Trim4_R2.fastq.gz -;
    """
}

process CUTADAPT_oligodT {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    label 'low'
    tag "${sample.sample}"

    input:
    val(index_step)
    tuple val (sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fastq.gz'), emit: trim_out
    
    script:
    """
    cutadapt -q 10 -a T{16} --minimum-length 46 --pair-filter=first --overlap 15 -o ${sample.sample}_R1_Ttrim.fastq.gz -p ${sample.sample}_Ttrim_R2.fastq.gz ${pe_reads[0]} ${pe_reads[1]} 
    """
}

//when read 1 has the cdna use after UMI_QSP
process CUTADAPT_QSP {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trim2.fastq.gz')
    script:
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[0]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample}_R2_trim2.fastq.gz -
    """
}

//command to trim the reads to a specific length
process CUTADAPT_short {

    container './containers/cutadapt.sif'
    publishDir "${params.outdir}/${index_step}_cutadapt_1", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trim2.fastq.gz')
    script:
    """
    cutadapt -u -28 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[0]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample}_R2_trim2.fastq.gz -
    """
}

//when read 2 has the cdna

process CUTADAPT_QSP_2 {

    container './containers/cutadapt.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_cutadapt_2", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trim2.fastq.gz')
    script:
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[1]} | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r2adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r2adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${sample}_R2_trim2.fastq.gz -
    """
}



// Trim Nextera adapters and Poly-A sequences from Read2 (the sequence of interest for alignment)

process CUTADAPT_QSP_2_v2 {

    container './containers/cutadapt.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_cutadapt_2", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trim2.fastq.gz')
    script:
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${pe_reads[1]} \
    | cutadapt -m 20 -O 3 --nextseq-trim=10 \
    -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" \
    | cutadapt -m 20 -O 3 -a "polyA=A{18}" \
    | cutadapt -m 20 -O 20 -g "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" \
    --discard-trimmed -o ${sample}_R2_trimmed.fastq.gz
    """
}

process CUTADAPT_targeted {

    container './containers/cutadapt.sif'
    label 'low'
    publishDir "${params.outdir}/${index_step}_cutadapt_2", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val(index_step)
    tuple val(sample), path(pe_reads)
    path primers_fasta

    output:
    tuple val(sample), file('*_trim2.fastq.gz')

    script:
    """
    cutadapt \
      -g file:${primers_fasta} \
      --max-n 0 \
      --discard-untrimmed \
      -o - ${pe_reads[1]} \
    | cutadapt \
        --nextseq-trim=10 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -a "A{20}" \
        -m 20 \
        -O 5 \
        -o ${sample}_R2_trim2.fastq.gz \
        -
    """
}



process TrimGalore {

    container 'quay.io/biocontainers/trim-galore:0.6.6--0'
    publishDir "${params.outdir}/${index_step}_trim_galore", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    val (index_step)
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*.fq.gz')

    script:
    """
    trim_galore --paired --quality 20 --length 20 --output_dir . ${pe_reads[0]} ${pe_reads[1]}
    mv ${sample}_R1_val_1.fq.gz ${sample}_trimmed_R1.fq.gz
    mv ${sample}_R2_val_2.fq.gz ${sample}_trimmed_R2.fq.gz
    """
}

