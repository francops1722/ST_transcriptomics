#!/bin/bash
#PBS -N TestElena_01
#PBS -l walltime=01:00:00
#PBS -l mem=64gb

cd /user/gent/446/vsc44685/ScratchVO_dir/ST_transcriptomics

ml Nextflow/23.04.2

export NXF_SINGULARITY_CACHEDIR="/user/gent/446/vsc44685/ScratchVO_dir/ST_transcriptomics_local"
export APPTAINER_TMPDIR="/user/gent/446/vsc44685/ScratchVO_dir/ST_transcriptomics_local"
export APPTAINER_CACHEDIR="/user/gent/446/vsc44685/ScratchVO_dir/ST_transcriptomics_local"

nextflow run NF_Elena.nf -c nextflow.config