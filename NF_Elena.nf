#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')

// input parameters 
params.reads ="/user/gent/446/vsc44685/DataVO_dir/Miseq_01122023/RawData/ERV_2/*_L001_{R1,R2}_001.fastq.gz"
// include processes
include {Subsample_Seqtk_or as subs} from "./modules/Seqtk" 
include {check_QC_or as check_QC_raw; check_QC_or as check_QC_raw_subs; check_QC_or as check_QC_trimmed} from "./modules/fastQC" 
include {check_QC_or as check_QC_trimmed2} from "./modules/fastQC" 
include {check_QC_or as check_QC_demux} from "./modules/fastQC" 
include {CUTADAPT_Trim1_or as trim} from "./modules/trimming"
include {UMI_extract as UMI_ext} from "./modules/UMItools"
include {Decode as demux} from "./modules/Demux"
include {CUTADAPT_Trim4_1 as trim2} from "./modules/trimming"
include {Star_Align_R2v2 as star} from "./modules/Star"
include {index_bam as index} from "./modules/Star"
include {index_bam_FC as index2} from "./modules/FeatureCounts"
include {sort_bam_FC as sort_bam} from "./modules/FeatureCounts"
include {UMI_count as count} from "./modules/UMItools"
include {FeatureCounts_BAM as FeatureCounts} from "./modules/FeatureCounts"
include {merge_Counts as merge} from "./modules/FeatureCounts"
include {merge_Counts as merge} from "./modules/FeatureCounts"
include {MULTIQC as check_star} from "./modules/fastQC"
include {MULTIQC as check_FC} from "./modules/fastQC"
include {plot_reads as plots} from "./modules/Make_plots"
include {Dedup_log as dedup_log} from "./modules/Log_dedup"


log.info """\
      LIST OF PARAMETERS
================================
            GENERAL
Data-folder      : $params.datadir
Results-folder   : $params.outdir
================================
          Subsampling
Number of reads  : $params.Subs
================================
          UMI-tools
UMI-type         : $params.UMItype          
================================          
      INPUT & REFERENCES 
Input-files      : $params.reads
Reference genome : $params.genome
GTF-file         : $params.gtf_file
================================
             STAR
Threads          : $params.threads
================================
"""

workflow {
    // set input data
    pe_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
    genome = file(params.genome)
    gtf = file(params.gtf_file)
    bar_file = file(params.bar_file)
    N = 3 //number of unique barcodes
    Nthr = 6 //6 when M=12 10 when M=36
    index_step = 0
    
    check_QC_raw("raw", pe_reads_ch)

    //is subsampling enabled?
    if (params.Subs) {
         subs(index_step, pe_reads_ch)
         pe_reads_ch = subs.out
         check_QC_raw_subs("raw_sub", pe_reads_ch)
         index_step += 1
    }     
    // run cutadapt
    trim(index_step, pe_reads_ch)
    //QC of trimmed reads
    check_QC_trimmed("trim", trim.out)
    //extract UMI
    index_step += 1
    UMI_ext(index_step, trim.out)
    //Demux barcodes
    index_step += 1
    demux(index_step, bar_file, N, Nthr, UMI_ext.out.fastq)
    check_QC_demux("demux", demux.out)
    //trim of poly A
    index_step += 1
    trim2(index_step, demux.out)
    //QC of trimmed reads
    check_QC_trimmed2("trim2", trim2.out)
    //STAR alignment of Read2
    index_step += 1
    star(index_step, trim2.out, genome)
    index(index_step, star.out.align_bam)
    //Feature Counts BAM input/output
    index_step += 1
    FeatureCounts(index_step, star.out.align_bam, gtf, index.out)
    //Sort BAM files output from FC
    sort_bam(index_step, FeatureCounts.out.FC_bam)
    index2(index_step, sort_bam.out.sorted_bam)
    FC_files = FeatureCounts.out.log_files.collect()
    check_FC("FCounts", FC_files)
    star_files = star.out.log_files.collect()
    check_star("star", star_files)
    //Counting molecules
    index_step += 1
    count(index_step, sort_bam.out.sorted_bam, index2.out)
    log_umi = count.out.log_files.collect()
    dedup_log(index_step, log_umi)
    //summary counts
    index_step += 1
    count_files = count.out.umi_counts.collect()
    merge(index_step, count_files)
    //make plots
    plots(params.outdir, merge.out.mock)
}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}