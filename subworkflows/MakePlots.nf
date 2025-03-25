#!/usr/bin/env nextflow
include {plot_reads; plot_map; plot_counts} from "../modules/Make_plots"

workflow All_Plots {
    take:
        files
        counts

    main:
        plot_reads(files, counts)
        plot_map(files, counts)
        plot_counts(files, counts)
    emit:
        plot_reads = plot_reads.out
        plot_map = plot_map.out
        plot_genecounts = plot_counts.out
}

