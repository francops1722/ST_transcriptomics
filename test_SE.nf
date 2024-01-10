#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')

// input parameters 
params.reads ="/user/gent/446/vsc44685/DataVO_dir/Miseq_01122023/RawData/02_27/*_L001_{R1,R2}_001.fastq.gz"
// include processes
include {Subsample_Seqtk_or as subs} from "./modules/Seqtk" 
include {UMI_QSP_2 as UMI_ext} from "./modules/UMItools" 
include {CUTADAPT_QSP_2 as trim} from "./modules/trimming"
include {Star_Align_QSP as star} from "./modules/Star"
include {index_bam as index} from "./modules/Star"
include {UMI_dedup_basic as dedup} from "./modules/UMItools"
include {FeatureCounts as FC} from "./modules/FeatureCounts"
include {merge_featureCounts as merge_FC} from "./modules/FeatureCounts"


workflow {
    // set input data
    pe_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
    genome = file(params.genome)
    gtf = file(params.gtf_file)
    //is subsampling enabled?
    if (params.Subs) {
         subs(pe_reads_ch)
         pe_reads_ch = subs.out
    } 
    UMI_ext(pe_reads_ch)
    trim(UMI_ext.out.fastq)  
    star(trim.out, genome) 
    index(star.out.align_bam) 
    dedup(star.out.align_bam, index.out)
    FC(dedup.out.dedup_files, gtf, index.out)
    input_files = FC.out.counts.collect()
    merge_FC(input_files)
}