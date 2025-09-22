#!/usr/bin/env nextflow
// include processes
include {Subsample_Seqtk_or as subs} from "../modules/Seqtk" 
include {UMI_QSP_2 as UMI_ext} from "../modules/UMItools" 
include {CUTADAPT_QSP_2 as trim} from "../modules/trimming"
include {CUTADAPT_targeted as trim_primer} from "../modules/trimming" //to trim primers from r2
include {Star_Align_QSP as star} from "../modules/Star" 
include {sort_bam_QSP as sort_bam} from "../modules/Star" 
// include {Star_Align_QSP_test as star} from "../modules/Star" //only for testing ffpe optimization
include {index_bam_QSP as index} from "../modules/Star"
// include {index_bam_QSP as index_trans} from "../modules/Star"
include {UMI_dedup_basic as dedup} from "../modules/UMItools"
// include {UMI_dedup_transcipt as dedup_trans} from "../modules/UMItools"
include {FeatureCounts as FC} from "../modules/FeatureCounts"
// include {Salmon_quant as salmon} from "../modules/UMItools"
include {merge_featureCounts as merge_FC} from "../modules/FeatureCounts"
// include {merge_TPMs as merge_tpm} from "../modules/FeatureCounts"

workflow QSPLike_workflow {
    // set input data
    pe_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
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

    if (params.primer_fasta) {
        primer_fasta = file(params.primer_fasta)
        trim_primer(index_step, UMI_ext.out.fastq, primer_fasta)
        trim_ch = trim_primer.out
    } else {
        trim(index_step, UMI_ext.out.fastq)
        trim_ch = trim.out
    }

    index_step += 1
    bam = star(index_step, trim_ch, genome)
    // sort_bam(index_step, star.out.align_bam)  
    bai = index(index_step, star.out.align_bam) 
    //bai_trans = index_trans(index_step, sort_bam.out.sorted_bam) 
    align_ch = star.out.align_bam.join(bai)
    // align_ch_trans = sort_bam.out.sorted_bam.join(bai_trans)
    index_step += 1
    dedup(index_step, align_ch)
    // dedup_trans(index_step, align_ch_trans)
    // // //dedup(index_step, align.out.bam_file)
    index_step += 1
    FC(index_step, dedup.out.dedup_files, gtf)
    // index_step += 1
    // salmon(index_step, dedup_trans.out.dedup_files, transcript_fasta)
    index_step += 1
    input_files = FC.out.counts.collect()
    // input_salmon_files = salmon.out.quant_folder.collect()
    // input_salmon_files.view()
    merge_FC(index_step, input_files)
    // merge_tpm(index_step, input_salmon_files)
}