#!/usr/bin/env nextflow
// include processes
include {Subsample_Seqtk_or as subs} from "../modules/Seqtk" 
include {UMI_QSP_2 as UMI_ext} from "../modules/UMItools" 
include {CUTADAPT_QSP_2 as trim} from "../modules/trimming"
include {Star_Align_QSP as star} from "../modules/Star"
include {index_bam_QSP as index} from "../modules/Star"
include {UMI_dedup_basic as dedup} from "../modules/UMItools"
include {FeatureCounts as FC} from "../modules/FeatureCounts"
include {merge_featureCounts as merge_FC} from "../modules/FeatureCounts"
include {Alignment_QSP as align} from "../subworkflows/StarWF"


workflow QSPLike_workflow {
    // set input data
    pe_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
    pe_reads_ch.view()
    // pe_reads_ch = Channel.fromFilePairs( params.reads , size:-1 , checkIfExists:true) { file -> file.name.replaceAll(".fastq.gz", '') }
    genome = file(params.genome)
    gtf = file(params.gtf_file)
    index_step = 0
    // is subsampling enabled?
    if (params.Subs) {
         subs(index_step, pe_reads_ch)
         pe_reads_ch = subs.out
         index_step += 1
    } 
    UMI_ext(index_step, pe_reads_ch)
    index_step += 1
    trim(index_step, UMI_ext.out.fastq)  
    // index_step += 1
    // bam = star(index_step, trim.out, genome)
    // // bam = star(index_step, pe_reads_ch, genome)  
    // bai = index(index_step, star.out.align_bam) 
    // align_ch = star.out.align_bam.join(bai).view()
    // //index_step += 1
    // dedup(index_step, align_ch)
    // // //dedup(index_step, align.out.bam_file)
    // index_step += 1
    // FC(index_step, dedup.out.dedup_files, gtf)
    // index_step += 1
    // input_files = FC.out.counts.collect()
    // merge_FC(index_step, input_files)
}