#!/usr/bin/env nextflow

// include processes
include {Subsample_Seqtk as subs} from "./modules/Seqtk" 
include {check_QC as QC_raw; check_QC as QC_raw_subs; check_QC as QC_trimmed} from "./modules/fastQC" 
include {check_QC_or as QC_trimmed2} from "./modules/fastQC" 
include {check_QC_or as check_QC_demux} from "./modules/fastQC" 
include {CUTADAPT_oligodT as trim} from "./modules/trimming"
include {UMI_extract as UMI_ext} from "./modules/UMItools"
include {Decode as demux} from "./modules/Demux"
include {CUTADAPT_Trim4 as trim2} from "./modules/trimming"
include {Star_Align_R2 as star} from "./modules/Star"
include {index_bam as index} from "./modules/Star"
include {index_bam_FC as index2} from "./modules/FeatureCounts"
include {sort_bam_FC as sort_bam} from "./modules/FeatureCounts"
include {FeatureCounts_BAM as FeatureCounts} from "./modules/FeatureCounts"
include {UMI_count as count} from "./modules/UMItools"
include {merge_Counts as merge} from "./modules/FeatureCounts"
include {MULTIQC as check_star} from "./modules/fastQC"
include {MULTIQC as check_FC} from "./modules/fastQC"


log.info """
=====================================================================
     _____ _______    ____                  _____  _   _          
    / ____|__   __|  / __ \\                |  __ \\| \\ | |   /\\    
   | (___    | |    | |  | |_ __   ___ ___ | |__) |  \\| |  /  \\   
    \\___ \\   | |    | |  | | '_ \\ / __/ _ \\|  _  /| . ` | / /\\ \\  
    ____) |  | |    | |__| | | | | (_| (_) | | \\ \\| |\\  |/ ____ \\ 
   |_____/   |_|     \\____/|_| |_|\\___\\___/|_|  \\_\\_| \\_/_/    \\_\\
=====================================================================  
                                                                                                                    
        LIST OF PARAMETERS
    ================================
                INPUT
    CSV-file         : $params.CSVfile 
    Results-folder   : $params.outdir
    ================================
            Subsampling
    Number of reads  : $params.Subs
    ================================
            UMI-tools
    UMI-type         : $params.UMItype          
    ================================   
            Decoding
    Length of Barcodes: $params.M
    ================================       
                STAR
    Threads          : $params.threads
    Reference genome : $params.genome
    GTF-file         : $params.gtf_file
    ================================
    """.stripIndent()

// Create the channels

samples_ch = Channel
                .fromPath(params.CSVfile)
                .splitCsv(header:true)

samples_ch.map{ 
    row -> id = row.Sample_id
    reads = tuple(row.R1_path, row.R2_path)
    N = row.N_barcodes
    return [[sample:id, type:N], reads]
    }.view()

bar_file = samples_ch.map{ row -> row.barcode_file} 
genome = file(params.genome)
gtf = file(params.gtf_file)

if (N == 42000) {
Nthr = params.N42
} else {
Nthr = params.N86
}


index_step = 0

workflow{
        
        QC_raw('raw', reads)

        if (params.Subs) {
                subs(index_step, sample_id, reads)
                reads_ch = subs.out
                QC_raw_subs('raw_subs', reads)
                index_step += 1
        } 

        // run cutadapt
        trim(index_step, sample_id, reads_ch)
        //extract UMI
        index_step += 1
        UMI_ext(sampleid, index_step, trim.out)
        //Demux barcodes
        index_step += 1
        demux(index_step, bar_file, N, Nthr, UMI_ext.out.fastq)
        check_QC_demux("demux", demux.out)
        //trim of poly A
        index_step += 1
        trim2(index_step, demux.out)
        QC_trimmed2('trim2', trim2.out)
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
        //QC reports
        FC_files = FeatureCounts.out.log_files.collect()
        check_FC("FCounts", FC_files)
        star_files = star.out.log_files.collect()
        check_star("star", star_files)
        //Counting molecules
        index_step += 1
        count(index_step, sort_bam.out.sorted_bam, index2.out)
        index_step += 1
        count_files = count.out.umi_counts.collect()
        merge(index_step, count_files)

}
        
workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



