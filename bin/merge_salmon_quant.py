#!/usr/bin/env python
import pandas as pd
import os
import sys

# Check if at least one folder is provided
if len(sys.argv) < 2:
    print("Usage: python merge_salmon_quant.py <folder1> [folder2 ...]")
    sys.exit(1)

# Get all folders from command line arguments
folders = sys.argv[1:]

# Initialize an empty dictionary to store dataframes
dfs = {}

# Process each folder
for folder in folders:
    quant_file = os.path.join(folder, 'quant.sf')
    if not os.path.exists(quant_file):
        print(f"quant.sf not found in {folder}")
        continue
        
    # Get the folder name (sample name)
    sample_name = os.path.basename(folder)
    
    # Read the file and select only Name and TPM columns
    df = pd.read_csv(quant_file, sep='\t', usecols=['Name', 'TPM'])
    
    # Rename TPM column to sample name
    df = df.rename(columns={'TPM': sample_name})
    
    # Store in dictionary
    dfs[sample_name] = df

if not dfs:
    print("No valid quant.sf files found in the provided folders")
    sys.exit(1)

# Merge all dataframes on Name column
merged_df = None
for sample_name, df in dfs.items():
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='Name', how='outer')

# Rename Name column to ensembl_transcript_id
merged_df = merged_df.rename(columns={'Name': 'ensembl_transcript_id'})

# Calculate sum of TPM values across all samples
merged_df['sum_TPM'] = merged_df.drop('ensembl_transcript_id', axis=1).sum(axis=1)

# Remove rows where sum_TPM = 0
merged_df = merged_df[merged_df['sum_TPM'] > 0]

# Remove the sum_TPM column
merged_df = merged_df.drop('sum_TPM', axis=1)

# Save to file
output_file = 'merged_salmon_TPM.txt'
merged_df.to_csv(output_file, sep='\t', index=False)
print(f"Results saved to {output_file}") 