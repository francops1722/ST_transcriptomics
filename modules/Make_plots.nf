#!/usr/bin/env nextflow


process plot_reads {
    publishDir "${params.outdir}/Plots", mode: 'copy', overwrite: true 
    label 'plots'

    input:
    path (files)
    path (counts)
    output:
    file('*')    
    
    script:
    """
    CountReads.r ${files}
    """
    }

process plot_map {
    publishDir "${params.outdir}/Plots", mode: 'copy', overwrite: true 
    label 'plots'

    input:
    path (files)
    path (counts)
    output:
    file('*')
    
    
    script:
    """
    MakeMap.r ${files} 
    """
    }

process plot_counts {
    publishDir "${params.outdir}/Plots", mode: 'copy', overwrite: true 
    label 'plots'

    input:
    path (files)
    path (counts)
    output:
    file('*')
    
    
    script:
    """
    MakeCountPlot.r ${files}
    """
    }



    process plot_dedup {
    publishDir "${params.outdir}/Plots", mode: 'copy', overwrite: true 
    label 'plots'

    input:
    path (files)
    path (counts)
    output:
    file('*')

    script:
    """
    DedupPlot.r ${files}
    """
    }

