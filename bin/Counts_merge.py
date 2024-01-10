#!/usr/bin/env python

import pandas as pd
import os
import re
import logging
import glob
import argparse


# Use re.search() to check if the string matches the pattern


def extract_sample_name(file_path):
    # Use regular expression to extract "Sample_SN" from the file path
    match = re.sub(r"_counts\.tsv\.gz$", "", file_path)
    if match:
        return match
    else:
        return "Unknown"

def merge_BarCounts(input_files, out_file):
    
    logger = logging.getLogger(__name__)
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    input_files = glob.glob("*_counts.tsv.gz")
    df_final = pd.DataFrame()
    i = 0

    for file in input_files:
        name = extract_sample_name(file)
        df_tmp = pd.read_table(file)
        df_tmp2 = df_tmp.groupby('cell').sum('count')
        df_tmp2 = df_tmp2.rename(columns={'count':name})
        if i == 0:
            df_final = df_tmp2.copy()
        else:
            df_final = df_final.merge(df_tmp2, on='cell', how = 'outer')
        i+=1
    
    logger.info("Writing to csv file {}".format(out_file))
    df_final = df_final.fillna(0)
    df_final.to_csv(out_file)

def merge_Counts(input_files, out_file):
    
    logger = logging.getLogger(__name__)
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    input_files = glob.glob("*_counts.tsv.gz")
    df_final = pd.DataFrame()
    i = 0

    for file in input_files:
        name = extract_sample_name(file)
        df_tmp = pd.read_table(file)
        df_tmp2 = df_tmp.groupby('gene').sum('count')
        df_tmp2 = df_tmp2.rename(columns={'count':name})
        if i == 0:
            df_final = df_tmp2.copy()
        else:
            df_final = df_final.merge(df_tmp2, on='gene', how = 'outer')
        i+=1
    
    logger.info("Writing to csv file {}".format(out_file))
    df_final = df_final.fillna(0)
    df_final.to_csv(out_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples in a project
    """)
    parser.add_argument("-Bar", "--results_file_name1",dest='out_file1', default='all_counts_bar.txt',
                                   help= "Name of the output barcode count file that will be created")
    parser.add_argument("-Count", "--results_file_name2",dest='out_file2', default='all_counts.txt',
                                   help= "Name of the output count file that will be created")
    parser.add_argument("-i", "--input_files", dest='input_files', metavar='<input_files>', nargs='+', default='*_counts.tsv.gz',
                                   help="Path to the inputfiles from FeatureCounts. ")
    args = parser.parse_args()

merge_BarCounts(args.input_files, args.out_file1)
merge_Counts(args.input_files, args.out_file2)
