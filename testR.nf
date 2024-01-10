#!/usr/bin/env nextflow

include {plot_reads as plots} from "./modules/Make_plots"


input_dir = "${params.outdir}/multiqc/**general_stats.txt"


workflow{
    files = Channel.fromPath(input_dir).collect()
    plots(files)
}
