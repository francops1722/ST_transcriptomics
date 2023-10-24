#!/usr/bin/env nextflow

// include processes
include {Subsample_Seqtk as subs} from "./modules/Seqtk" 
include {check_QC as QC_raw; check_QC as QC_trimmed} from "./modules/fastQC" 
include {check_QC as QC_trimmed2} from "./modules/fastQC" 
include {CUTADAPT_oligodT as trim} from "./modules/trimming"
include {UMI_extract as UMI_ext} from "./modules/UMItools"
include {Decode as demux} from "./modules/Demux"
include {CUTADAPT_Trim4 as trim2} from "./modules/trimming"
include {Star_Align_R2 as star} from "./modules/Star"
include {index_bam as index} from "./modules/Star"
include {index_bam_FC as index2} from "./modules/FeatureCounts"
include {sort_bam_FC as sort_bam} from "./modules/FeatureCounts"
include {FeatureCounts_BAM as FeatureCounts} from "./modules/FeatureCounts"
include {UMI_count as count} from "./modules/FeatureCounts"


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

sample_id = samples_ch.map{ row -> row.Sample_id} 
reads_ch = samples_ch.map{ row -> tuple(row.R1_path, row.R2_path) }
N = samples_ch.map{ row -> row.N_barcodes} 
bar_file = samples_ch.map{ row -> row.barcode_file} 
genome = file(params.genome)
gtf = file(params.gtf_file)

if (N == 42000) {
Nthr = params.N42
} else {
Nthr = params.N86
}

workflow {
        if (params.Subs) {
        subs(sample_id, reads_ch)
        reads_ch = subs.out
        } 
        QC_raw('raw', reads_ch)
        // run cutadapt
        trim(sample_id, reads_ch)
        //trim of poly A
        trim2(trim.out)
        QC_trimmed2('trim', trim2.out)
        //extract UMI
        UMI_ext(trim2.out)
        //Demux barcodes
        demux(bar_file, N, Nthr, UMI_ext.out.fastq)
        //STAR alignment of Read2
        star(demux.out, genome)
        index(star.out.align_bam)
        //Feature Counts BAM input/output
        FeatureCounts(star.out.align_bam, gtf, index.out)
        //Sort BAM files output from FC
        sort_bam(FeatureCounts.out.FC_bam)
        index2(sort_bam.out.sorted_bam)
        //Counting molecules
        count(sort_bam.out.sorted_bam, index2.out)
        
}
workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



