#!/usr/bin/env nextflow
include { validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-validation'
include {check_QC as QC_raw; check_QC as QC_raw_subs; check_QC as QC_trim2;check_QC as check_QC_demux} from "../subworkflows/CheckQC"
include {Subsample_Seqtk as subs} from "../modules/Seqtk" 
include {CUTADAPT_oligodT as trim; CUTADAPT_Trim2 as trim2} from "../modules/trimming"
include {UMI_extract as UMI_ext} from "../modules/UMItools"
include {Decode_batch_4 as demux_batch} from "../modules/barcode_demux/main"
include {Star_Align_R2v2 as star} from "../modules/Star"
include {index_bam as index} from "../modules/Star"
include {FeatureCounts_BAM as FeatureCounts} from "../modules/FeatureCounts"
include {index_bam_FC as index2} from "../modules/FeatureCounts"
include {sort_bam_FC as sort_bam} from "../modules/FeatureCounts"
include {MULTIQC as check_star} from "../modules/fastQC"
include {MULTIQC as check_FC} from "../modules/fastQC"
include {UMI_count as count} from "../modules/UMItools"
include {merge_Counts as merge; table_reads as reads} from "../modules/FeatureCounts"
include {All_Plots as plots} from "../subworkflows/MakePlots"

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
    Reference genome : $params.genome
    GTF-file         : $params.gtf_file
    ================================
    """.stripIndent()

// Create the channels
samples_ch = Channel
                .fromPath(params.CSVfile)
                .splitCsv(header:true)

data = samples_ch.map{ 
    row -> id = row.Sample_id
    reads = tuple(row.R1_path, row.R2_path)
    N = row.N_barcodes
    barcode_file = row.barcode_file
    return [[sample:id, type:N, bar:barcode_file], reads]
    }



genome = file(params.genome)
gtf = file(params.gtf_file)

workflow ST_workflow {
    index_step = 0
    QC_raw('raw', data)
    
    if (params.Subs) {
        subs(index_step, data)
        data = subs.out
        QC_raw_subs('raw_subs', data)     
        } 
    
    // run cutadapt
    index_step += 1
    trim(index_step, data)
    index_step += 1 
    UMI_ext(index_step, trim.out)
    //Demux barcodes
    index_step += 1
    demux_batch(index_step, UMI_ext.out.fastq)
    check_QC_demux("demux", demux_batch.out) 
    // trim of poly A
    index_step += 1
    trim2(index_step, demux_batch.out)
    QC_trim2('trim2', trim2.out)
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
    reads(index_step, params.outdir, merge.out.mock)
    //make plots
    plots(params.outdir, reads.out)
}

