#!/bin/bash
# avoid that Apptainer uses $HOME/.cache
#it has to be set to a dir with enough space for all the image asociated files
export APPTAINER_CACHEDIR=/tmp/$USER/apptainer/cache 
# instruct Apptainer to use temp dir on local filessytem
export APPTAINER_TMPDIR=/tmp/$USER/apptainer/tmpdir
# specified temp dir must exist, so create it
mkdir -p $APPTAINER_TMPDIR

# Define the target folder for Singularity images
output_folder="./containers"
mkdir -p $output_folder
# List of Docker images to download
docker_images=(
    "docker://quay.io/biocontainers/cutadapt:4.5--py39hf95cd2a_0"
    "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    "docker://ewels/multiqc:latest"
    "docker://quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    "docker://quay.io/staphb/samtools:1.16"
    "docker://quay.io/biocontainers/umi_tools:1.1.4--py310h4b81fae_2"
    "docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    "docker://quay.io/biocontainers/star:2.7.8a--0"
)
# Loop through the list and convert each Docker image to Singularity

for image in "${docker_images[@]}"; do
  # Extract the image name (last part) to use as the Singularity image name
  image_name=$(echo "$image" | awk -F'[:/]' '{print $(NF-1)}')
  # Pull the Docker image and convert it to Singularity
  echo "Building $image_name"
  apptainer pull "$output_folder/$image_name.sif" "$image"
done

echo "Singularity images downloaded and saved to $output_folder"
