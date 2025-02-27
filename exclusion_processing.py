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

#print("\nMerged with specific columns:")
#print(merged_data_specific.head())
#################
# code for is_effective == t and exclusion_ids 1,7,9,14,29,94,96,99

# Usage
ids_to_keep = [29, 94, 96, 99]
#ids_to_omit = [1, 2, 7, 9, 14, 26] #26514
ids_to_omit = [1, 7, 9, 14, 2] #26546

assembly_level_to_exclude = ['partial']
excluded_protein_counts = [0]

###############################################################################
def process_exclusion_data(df, 
                           ids_to_keep, 
                           ids_to_omit, 
                           excluded_protein_counts=[0], 
                           assembly_level_to_exclude = ['partial', 'full']
                           ):
    """
    Processes exclusion reasons from a dataset, handling cases where multiple exclusion reasons
    may exist for the same UPID, applying logic to keep or omit records based on user-defined criteria.

    Args:
    - file_path (str): Path to the TSV file containing the data.
    - ids_to_keep (list): List of exclusion IDs that, if found, will lead to keeping the record.
    - ids_to_omit (list): List of exclusion IDs that, if found, will lead to omitting the record.
    - excluded_protein_counts (list): List of protein counts that should lead to complete removal of the record.
    - check_immunity (bool): Whether to check the 'exclusion_immunity' field before processing.

    Returns:
    - pd.DataFrame: Processed DataFrame with UPIDs and their exclusion reasons collapsed, filtered by the provided criteria.
    - pd.DataFrame: Log DataFrame containing details of any special conditions or flags.
    """
    # Load the data file
    #df = pd.read_csv(file_path, sep='\t')

    # Initialize log data structure
    #log_data = {
    #    'multiple_ids_to_keep': [],
    #    'multiple_ids_to_omit': [],
    #    'conflicted_ids': []
    #}


    print(f"\nExcluding proteomes with protein count: {excluded_protein_counts}")
    df2 = df[~df['protein_count'].isin(excluded_protein_counts)]
    print(len(df2))
    
    print(f"\nExcluding proteomes with assembly level: {assembly_level_to_exclude}")
    df2 = df2[~df2['assembly_level'].isin(assembly_level_to_exclude)]
    print(len(df2))
    
    
    
############
# Custom aggregation function to concatenate values
def concatenate_values(series):
    return ', '.join(map(str, series.unique()))  # Concatenate unique values as strings

# Custom aggregation function for 'set'
def to_set(series):
    return set(series)

# reorder cols
def move_cols_to_end(df, cols):
    return df[[x for x in df.columns if not x in cols] + cols]

############
df = merged_data_all.copy()
print(f"\nLength of input data: {len(df)}")
print(f"\nCount effective exclusions:\n {df['is_effective'].value_counts()}")

print(f"\nExcluding proteomes with protein count: {excluded_protein_counts}")
df2 = df[~df['protein_count'].isin(excluded_protein_counts)]
print(f"\nLength of df: {len(df2)}")

print(f"\nExcluding proteomes with assembly level: {assembly_level_to_exclude}")
df2 = df2[~df2['assembly_level'].isin(assembly_level_to_exclude)]
print(f"\nLength of df: {len(df2)}")

#3 upids belongig to 10 records have partial
#	upid
#3944	UP000243639
#3945	UP000243639
#3946	UP000243639
#4450	UP000243565
#4451	UP000243565
#4452	UP000243565
#26283	UP000242121
#26284	UP000242121
#26285	UP000242121
#26304	UP000464300

print(f"\nCount of effective exclusions with 'protein_count > 0' and 'assembly level not partial':\n {df2['is_effective'].value_counts()}")


# Step 1: Identify upids to remove
grouped = df2.groupby('upid')['exclusion_id'].agg(set)
#grouped2 = df2.groupby('upid')['exclusion_id'].agg(list)

upids_to_remove = grouped[
    grouped.apply(lambda ids: any(i in ids_to_omit for i in ids) and not any(i in ids_to_keep for i in ids))
].index

upids_to_filter = grouped[
    grouped.apply(lambda ids: any(i in ids_to_keep for i in ids) and any(i in ids_to_omit for i in ids))
].index

a = upids_to_remove.to_list() + upids_to_filter.to_list()

check2 =df2[~df2['upid'].isin(upids_to_remove.to_list() + upids_to_filter.to_list())]
#some upids to check: UP000000685, UP000002642, UP000235432 (8772)

#
# Step 2: Filter DataFrame
df_filtered = df2[~df2['upid'].isin(upids_to_remove)]  # Remove unwanted upids
df_filtered = df_filtered[
    df_filtered['upid'].isin(upids_to_filter) & df_filtered['exclusion_id'].isin(ids_to_keep) | 
    ~df_filtered['upid'].isin(upids_to_filter)  # Keep only keep IDs where needed
]

# Step 3: Define aggregation functions for each column
agg_funcs = {
    'exclusion_id': to_set,  # Use a custom function (to_set) for 'exclusion_id'
    'id_description': concatenate_values,  # Concatenate unique values for id_description
    'is_effective': concatenate_values,  # Concatenate unique values for is_effective
    'source': concatenate_values,  # Concatenate unique values for source
}

# For all other columns, use 'first' as the aggregation function
other_cols = df2.columns.difference(['upid', 'protein_count', 'exclusion_id']).tolist()

# Default aggregation function for other columns
for col in other_cols:
    if col not in agg_funcs:  # Ensure we don't overwrite custom columns
        agg_funcs[col] = 'first'

# Step 4: Apply groupby with custom aggregation functions
df_grouped = df_filtered.groupby(['upid', 'protein_count']).agg(agg_funcs).reset_index()

# Output result
print(df_grouped)

# To check for duplicates in `df_grouped`
df_grouped['protein_count'].describe()
df_grouped['upid'].duplicated().value_counts()

df_filtered['upid'].nunique()
df_grouped['upid'].nunique()

# Check for any remaining duplicates after grouping
duplicates = df_grouped[df_grouped.duplicated(subset=['upid', 'protein_count'], keep=False)]
print(duplicates)

df2['is_effective'].value_counts()
df_grouped['is_effective'].value_counts()
df_grouped['protein_count'].describe()
######

cols_move_end = ['exclusion_id','id_description','is_effective','source']
df_grouped_reordered = move_cols_to_end(df_grouped, cols_move_end)


outfile_df = basedir + "/spneumo_biosample_info_processed.out"
df_grouped_reordered.to_csv(outfile_df, sep="\t", index=False )

outfile_ids = basedir + "/proc-ids.txt"
df_grouped_reordered['upid'].to_csv(outfile_ids, index=False)

#############################

