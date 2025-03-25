#!/bin/bash

# Check if directory is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory> <output>"
    exit 1
fi

# Directory containing files
directory="$1"
output_file="$1/summary.csv"

# Check if directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# Output file
#output_file="/user/gent/446/vsc44685/ScratchVO_dir/ST_800K/Dedup_basic_extra/SummaryDedup.csv"

# Remove output file if exists
rm -f "$output_file"

# Loop through each file in the directory
for file in "$directory"/*.log; do
    if [ -f "$file" ]; then
        # Extract values from file
        input_reads=$(grep -oP 'Input Reads: \K.*' "$file")
        post_dedup_st=$(grep -oP 'Number of \(post deduplication\) reads counted: \K.*' "$file")
        post_dedup_qsp=$(grep -oP 'Number of reads out: \K.*' "$file")
        
        # Append filename and values to output file
        echo "$(basename "$file"), $input_reads, $post_dedup_st, $post_dedup_qsp" >> "$output_file"
    fi
done

echo "Extraction complete. Output saved to '$output_file'."
