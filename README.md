# Spatial Transcriptomics OncoRNA Lab
pipeline for the analysis of Spatial Transcriptomics data generated from our custom microarray

## Requirements
1. singularity/apptainer
2. Nextflow
3. Pytorch
4. Nvidia Cuda

## Instalation

1. clone the repository
`git clone https://github.com/francops1722/ST_transcriptomics.git`
2. Download containers
`./build_containers.sh`
4. Add all need paramters in the csvParams.config
5. Load nextflow module
`ml Nextflow/24.10.2`
7. Run pipeline
`nextflow run ST_analysis.nf -c nextflow.config -with-dag ST_flowchart.png`

