#!/usr/bin/env nextflow

include {plot_reads as plots} from "./modules/Make_plots"

workflow{
    plots(params.outdir, "test")
}
