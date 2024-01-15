#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')

// input parameters 
params.reads ="/user/gent/446/vsc44685/DataVO_dir/QSP_MAQCA/Raw/*_{R1,R2}.fastq.gz"
// include processes
include {Subsample_Seqtk_or as subs} from "./modules/Seqtk" 
include {UMI_QSP_2 as UMI_ext} from "./modules/UMItools" 
include {CUTADAPT_QSP as trim} from "./modules/trimming"
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
    index_step = 0
    //is subsampling enabled?
    if (params.Subs) {
         subs(index_step, pe_reads_ch)
         pe_reads_ch = subs.out
         index_step += 1
    } 
    UMI_ext(index_step, pe_reads_ch)
    index_step += 1
    trim(index_step, UMI_ext.out.fastq)  
    index_step += 1
    star(index_step, trim.out, genome) 
    index(index_step, star.out.align_bam) 
    index_step += 1
    dedup(index_step, star.out.align_bam, index.out)
    index_step += 1
    FC(index_step, dedup.out.dedup_files, gtf, index.out)
    index_step += 1
    input_files = FC.out.counts.collect()
    merge_FC(index_step, input_files)
}