#!/bin/bash

# Check if directory is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Directory containing BAM files
directory="$1"

# Check if directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# Create logs directory if it doesn't exist
logs_directory="logs_dedup"
mkdir -p "$logs_directory"

# Loop through each BAM file in the directory
for bam_file in "$directory"/*sorted.bam; do
    if [ -f "$bam_file" ]; then
        # Get sample name from BAM file
        sample=$(basename "$bam_file" .bam)

        # Run umi_tools dedup command
        umi_tools dedup -I "$bam_file" -S "${sample}_dedup.bam" \
            --multimapping-detection-method=NH \
            --output-stats="${logs_directory}/${sample}_deduplicated.txt" \
            --log="${logs_directory}/${sample}_deduplication.log"

        echo "Deduplication complete for $sample."
    fi
done

echo "All deduplications complete."
