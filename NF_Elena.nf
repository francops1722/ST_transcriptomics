#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')

// input parameters 
params.reads ="/user/gent/446/vsc44685/DataVO_dir/Miseq_01122023/RawData/ERV_2/*_L001_{R1,R2}_001.fastq.gz"
// include processes
include {Subsample_Seqtk_or as subs} from "./modules/Seqtk" 
include {check_QC_or as check_QC_raw; check_QC_or as check_QC_trimmed} from "./modules/fastQC" 
include {check_QC_or as check_QC_trimmed2} from "./modules/fastQC" 
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
    //is subsampling enabled?
    if (params.Subs) {
         subs(pe_reads_ch)
         pe_reads_ch = subs.out
    }     
    
    //pass the 'step' and the raw subs reads to the QC subworkflow
    check_QC_raw("raw_sub", pe_reads_ch)
    // run cutadapt
    trim(pe_reads_ch)
    //QC of trimmed reads
    check_QC_trimmed("trim", trim.out)
    //extract UMI
    UMI_ext(trim.out)
    //Demux barcodes
    demux(bar_file, N, Nthr, UMI_ext.out.fastq)
    //trim of poly A
    trim2(demux.out)
    //QC of trimmed reads
    check_QC_trimmed2("trim2", trim2.out)
    //STAR alignment of Read2
    star(trim2.out, genome)
    index(star.out.align_bam)
    //Feature Counts BAM input/output
    FeatureCounts(star.out.align_bam, gtf, index.out)
    //Sort BAM files output from FC
    sort_bam(FeatureCounts.out.FC_bam)
    index2(sort_bam.out.sorted_bam)
    //Counting molecules
    count(sort_bam.out.sorted_bam, index2.out)
    input_files = count.out.collect()
    merge(input_files)
}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}