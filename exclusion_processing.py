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
basedir =  homedir + "/Documents/arise/spneumo_dataset"
#basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

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
###############################################################################
import pandas as pd
from tqdm.auto import tqdm

df = merged_data_all.copy()

def process_exclusion_data(df, ids_to_keep, ids_to_omit, excluded_protein_counts, check_immunity=True):
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
    log_data = {
        'multiple_ids_to_keep': [],
        'multiple_ids_to_omit': [],
        'conflicted_ids': []
    }

    # Exclude records by protein count
   df2 = df[~df['protein_count'].isin(excluded_protein_counts)]

    # Group by 'upid' and apply processing
#    tqdm.pandas(desc="Processing groups")

    # return cleaned_df, log_df
    #df.groupby(['is_excluded']).size()
    df2.groupby(['exclusion_immunity']).size()
    df2.groupby(['exclusion_immunity', 'is_excluded']).size()
    df2.groupby(['is_excluded', 'exclusion_immunity']).size()
    df2.groupby(['is_excluded', 'is_effective']).size()


    #condition1 = df2['exclusion_id'] == 't'
    #condition2 = df2['exclusion_immunity'] == 't'
    condition_keep = df2['exclusion_id'].isin(ids_to_keep) 
    condition_omit = df2['exclusion_id'].isin(ids_to_omit) 
    
        
    #df3 = df2[condition1 | condition2] 
    #df3 = df2[condition1 & condition2 & condition3] 
    len(df2[condition_omit])
    df3 = df2[~condition_omit]
    len(df3) == len(df2) - len(df2[condition_omit])
    
    df3.groupby('upid').size()
    
    # count duplicates of upid
    dups = df3[df3.duplicated(subset='upid')]
    
    df3[condition_keep]
    df4 = df3[condition_keep]
    
    dups = df4[df4.duplicated(subset='upid')]

    a = df[df['upid'] == 'UP000304910']
    b = df[df['upid'].isin(list(dups['upid']))]
    c = b.groupby('upid').size()
    df4.groupby(['is_excluded', 'is_effective']).size()
    
    
    # Select all columns except 'exclusion_id' for duplicate checking
    cols_to_checkL = df.columns.difference(['exclusion_id', 'biosample']).tolist()


#########
    # Identify duplicates based on 'upid' and all other columns except 'exclusion_id'
    df.duplicated(subset=['upid', 'protein_count']).value_counts()
    
    df.duplicated(subset=cols_to_checkL).value_counts()
    
    df2.duplicated(subset=['upid', 'protein_count'], keep=False).value_counts()
    
    df2['is_duplicate'] = df2.duplicated(subset=['upid', 'protein_count'], keep=False)
    df2['is_duplicate'].value_counts()

    
    #if df2['is_duplicate'].any():
    #    print(df2['is_duplicate'].value_counts())
        
    
###########    

    df2['is_duplicate'] = df2.duplicated(subset=['upid', 'protein_count'], keep=False)
    df2['is_duplicate'].value_counts()
    
    # Create the condition to identify rows to omit
    condition_omit = (df2['is_duplicate'] == True) & (df2['exclusion_id'].isin(ids_to_omit))
    
    # Use the condition to subset the DataFrame by keeping rows that do not meet the condition
    filtered_df1 = df2[~condition_omit]
    
    print(filtered_df1)


    # Calculate total rows in the original DataFrame
    total_rows = len(df2)
    
    # Create the condition to identify rows to omit
    condition_omit = (df2['is_duplicate'] == True) & (df2['exclusion_id'].isin(ids_to_omit))
    
    # Count the number of rows that meet the condition to omit
    rows_to_omit = df2[condition_omit].shape[0]
    
    # Calculate the expected length of the final DataFrame
    expected_length_filtered_df1 = total_rows - rows_to_omit
    
    # Apply the condition to subset the DataFrame
    filtered_df1 = df2[~condition_omit]
    
    # Actual length of the final DataFrame
    actual_length_filtered_df1 = len(filtered_df1)
    
    # Output the expected and actual lengths for verification
    print("Expected length of final DataFrame:", expected_length_filtered_df1)
    print("Actual length of final DataFrame:", actual_length_filtered_df1)
    
    try:
        # Assert to check if the actual length matches the expected length
        assert actual_length_filtered_df1 == expected_length_filtered_df1, \
            "Mismatch in expected and actual lengths of first filtered DataFrame"
        print("Success: The lengths match. Filtering process is verified.")
    except AssertionError as e:
        print(e)

###
    
    # Assuming 'filtered_df1' is a predefined DataFrame from prior steps
    
    # Calculate total rows in the original DataFrame
    total_rows = len(filtered_df1)
    
    # Create the condition to identify rows to keep (not omit)
    #condition_keep = (filtered_df1['is_duplicate'] == True) & (filtered_df1['exclusion_id'].isin(ids_to_keep))
    condition_keep = (filtered_df1['exclusion_id'].isin(ids_to_keep))

    # Count the number of rows that meet the condition to keep
    rows_to_keep = filtered_df1[condition_keep].shape[0]
    
    # Calculate the expected length of the final DataFrame after keeping
    expected_length_filtered_df2 = rows_to_keep
    
    # Apply the condition to subset the DataFrame further if needed
    # Assuming here you want to keep these rows:
    filtered_df2  = filtered_df1[condition_keep]
    
    # Actual length of the final DataFrame
    actual_length_filtered_df2 = len(filtered_df2)
    
    # Output the expected and actual lengths for verification
    print("Expected length of final DataFrame:", expected_length_filtered_df2)
    print("Actual length of final DataFrame:", actual_length_filtered_df2)
    
    # Final Sanity Check
    #assert actual_length_filtered_df2 == expected_length_filtered_df2, "Mismatch in expected and actual lengths of final DataFrame"

    # Sanity Check
    try:
        # Assert to check if the actual length matches the expected length
        assert actual_length_filtered_df2 == expected_length_filtered_df2, \
            "Mismatch in expected and actual lengths of first filtered DataFrame"
        print("Success: The lengths match. Filtering process is verified.")
    except AssertionError as e:
        print(e)

##
filtered_df2['is_duplicate'].value_counts()
filtered_df2['upid'].nunique()

df4[df4.duplicated(subset= ['upid', 'protein_count'])]
df4[df4.duplicated(subset= ['upid'])]

df4['is_duplicate'].value_counts()
df4['upid'].nunique()

###

filtered_df2['upid'].nunique()
filtered_df2['exclusion_id'].isin(ids_to_keep).value_counts()
filtered_df2['exclusion_id'].isin(ids_to_omit).value_counts()

filtered_df2.groupby(['exclusion_id']).size()
filtered_df2['upid'].nunique() #17911


filtered_df2 = filtered_df2[filtered_df2['assembly_level'] !='partial']


###
dups_df2 = filtered_df2[filtered_df2.duplicated(subset= ['upid', 'protein_count'])]
check = filtered_df2[filtered_df2['upid'].isin(dups_df2['upid'])]
check['protein_count'].value_counts()
###

# Define custom aggregation for specific columns
def concatenate_values(series):
    return ','.join(map(str, series))  # Join all values into a comma-separated string without set conversion

# Set up aggregation rules
aggregations = {
    'exclusion_id': concatenate_values,
    'id_description': concatenate_values,
    'is_effective': concatenate_values,
    'source': concatenate_values,
}

# Default all other columns to use 'first' if not specifically mentioned
for col in filtered_df2.columns:
    if col not in aggregations:
        aggregations[col] = 'first'

# Group by 'upid' and aggregate
aggregated_df = filtered_df2.groupby('upid').agg(aggregations)

aggregated_df = aggregated_df.reset_index(drop=True)

#aggregated_df = aggregated_df.reset_index()  # 'upid' becomes a column here properly
# Display the result
print(aggregated_df)
a = aggregated_df[aggregated_df['upid'] == 'UP000304910']

# reorder cols
def move_cols_to_end(df, cols):
    return df[[x for x in df.columns if not x in cols] + cols]

cols_move_end = ['exclusion_id','id_description','is_effective','source']

aggregated_df_reordered = move_cols_to_end(aggregated_df, cols_move_end)

outdf = aggregated_df_reordered.drop(columns=['is_duplicate'])


outfile_df = basedir + "/spneumo_biosample_info_processed.out"
outdf.to_csv(outfile_df, sep="\t", index=False )

outfile_ids = basedir + "/proc-ids.txt"
outdf['upid'].to_csv(outfile_ids, index=False)

outdf['protein_count'].describe()

#############################


# Usage
ids_to_keep = [29, 94, 96, 99]
ids_to_omit = [1, 2, 7, 9, 14]
assembly_levels_to_exclude = ['partial']
excluded_protein_counts = [0]
cleaned_data, logs = process_exclusion_data(merged_data_all
                                            , ids_to_keep
                                            , ids_to_omit
                                            , excluded_protein_counts
                                            , assembly_levels_to_exclude)

print(cleaned_data.head())
print(logs)

