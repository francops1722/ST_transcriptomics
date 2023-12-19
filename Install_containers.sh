#!/bin/bash
# Create the output folder if it doesn't exist
mkdir -p "./containers"
# Define the target folder for Singularity images
output_folder="./containers"
# List of Docker images to download
docker_images=(
    "docker://quay.io/biocontainers/cutadapt:4.5--py39hf95cd2a_0"
    "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    "docker://ewels/multiqc:latest"
    "docker://quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    "docker://quay.io/biocontainers/samtools:1.2--0"
    "docker://quay.io/biocontainers/umi_tools:1.1.4--py310h4b81fae_2"
    "docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    "docker://quay.io/biocontainers/star:2.7.8a--0"
)

# Loop through the list and convert each Docker image to Singularity
for image in "${docker_images[@]}"; do
  # Extract the image name (last part) to use as the Singularity image name
  image_name=$(echo "$image" | awk -F'/' '{print $NF}')
  # Pull the Docker image and convert it to Singularity
  singularity pull "$output_folder/$image_name.sif" "$image"
done

echo "Singularity images downloaded and saved to $output_folder"
