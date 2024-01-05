#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')

// input parameters 
params.reads ="/user/gent/446/vsc44685/ScratchVO_dir/Out_test3/Counts/*"

include {merge_Counts as merge} from "./modules/FeatureCounts"

workflow {
    files = Channel.fromPath(params.reads)
    list = files.collect()
    merge(list)
}
