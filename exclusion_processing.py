#!/bin/env python
import sys, os
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib_venn import venn2
#from matplotlib.patches import Patch
#import seaborn as sns
#import random
#import argparse
#from tqdm import tqdm

######
# EXCLUSION processing
#if a up proteome is excluded (is_excluded == t), then look for exclusion_immunity. :
#    If exclusion_immunity==t, keep it.
#    If exclusion_immunity==f, move on to exclusion_id.
#    If exclusion_id = 2, then omit # partial
#    If exclusion_id = 1, then omit # genome length too large
#    If exclusion_id = 7, then omit # genome length too small
######
# Read data
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
up_spneumo_proteomes = basedir + "/spneumo_biosample_info.out"
exclusion_reasons = basedir + "/exclusion_reasons.tsv"

# Read the data
#df_up_original = pd.read_csv(up_file, sep="\t")
#print(f"\nDF shape: {df_up_original.shape}")
#print(f"\nLength of UP DF: {len(df_up_original)}")

# Step 0: omit the ones with protein_count == 0
#df_up_all = df_up_original[df_up_original['protein_count']> 0]
#print(f"\nNo. of of entries with 0 proteins: {len(df_up_original[df_up_original['protein_count']== 0])}")
#print(f"\nLength of UP DF after removing proteomes with 0 proteins: {len(df_up_all)}")
   
###############################################################################
def merge_tsv_files(file1_path, file2_path, 
                    file1_merge_col='exclusion_id', 
                    file2_merge_col='id', 
                    include_columns="All",
                    join_type='left'):
    """
    Merges two TSV files based on a merge column. Allows for different column names
    in each file for the merge key and includes specified columns from the second file or all by default.
    The type of join performed can be specified to control how rows from each dataset are included in the result.

    Args:
    - file1_path (str): Path to the first TSV file.
    - file2_path (str): Path to the second TSV file, containing additional data like exclusion meanings.
    - file1_merge_col (str): Column name in the first file to merge on. Default is 'exclusion_id'.
    - file2_merge_col (str): Column name in the second file to merge on. Default is 'id'.
    - include_columns (str or list): Specific columns from the second file to add to the first. If "All", all columns are added.
    - join_type (str): Type of join to perform during the merge. Options include 'left', 'right', 'outer', 'inner'. Default is 'left'.

    Returns:
    - pd.DataFrame: A DataFrame with the merged content.

    Raises:
    - ValueError: If any specified columns in `include_columns` do not exist in the second file.
    """
    
    # Load the data files
    df1 = pd.read_csv(file1_path, sep='\t')
    df2 = pd.read_csv(file2_path, sep='\t')

    # Rename the merge column in file2 if necessary
    if file2_merge_col != file1_merge_col:
        df2.rename(columns={file2_merge_col: file1_merge_col}, inplace=True)

    # Check if specific columns are to be included from the second file
    if include_columns != "All":
        if not set(include_columns).issubset(df2.columns):
            missing_cols = set(include_columns) - set(df2.columns)
            raise ValueError(f"The following columns are specified in include_columns but do not exist in second file {missing_cols}")
        columns_to_include = [file1_merge_col] + include_columns
        df2 = df2[columns_to_include]

    # Merge the dataframes on the specified column
    merged_df = pd.merge(df1, df2, on=file1_merge_col, how=join_type)

    return merged_df
###############################################################################
# Example usage:
# Merging and including all columns from the second file (default behavior)

merged_data_all = merge_tsv_files(file1_path=up_spneumo_proteomes, 
                                  file2_path=exclusion_reasons,
                                  include_columns="All")

# Merging and specifying only certain columns to include from the second file
merged_data_specific = merge_tsv_files(file1_path=up_spneumo_proteomes, 
                                  file2_path=exclusion_reasons,
                                  include_columns=['id_description'])

print("Merged with all columns (default):")
print(merged_data_all.head())

print("\nMerged with specific columns:")
print(merged_data_specific.head())
