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
import pandas as pd

def merge_data(file1, file2, 
                    file1_merge_col='exclusion_id', 
                    file2_merge_col='id', 
                    include_columns="All",
                    join_type='left',
                    delimiter = "\t",
                    drop_cols=[]):
    """
    Merges two datasets based on a merge column. Allows for different column names
    in each dataset for the merge key and includes specified columns from the second dataset or all by default.
    The type of join performed can be specified to control how rows from each dataset are included in the result.
    Each dataset can either be a DataFrame or a file path to a TSV file.

    Args:
    - file1 (str or pd.DataFrame): First dataset or path to the TSV file.
    - file2 (str or pd.DataFrame): Second dataset or path to the TSV file, containing additional data like exclusion meanings.
    - file1_merge_col (str): Column name in the first dataset to merge on. Default is 'exclusion_id'.
    - file2_merge_col (str): Column name in the second dataset to merge on. Default is 'id'.
    - include_columns (str or list): Specific columns from the second dataset to add to the first. If "All", all columns are added.
    - join_type (str): Type of join to perform during the merge. Options include 'left', 'right', 'outer', 'inner'. Default is 'left'.

    Returns:
    - pd.DataFrame: A DataFrame with the merged content.

    Raises:
    - ValueError: If any specified columns in `include_columns` do not exist in the second dataset.
    """

    # Load or use data frames
    if isinstance(file1, str):
        df1 = pd.read_csv(file1, sep=delimiter)
    elif isinstance(file1, pd.DataFrame):
        df1 = file1
    else:
        raise ValueError("file1 must be a filepath or a DataFrame")

    if isinstance(file2, str):
        df2 = pd.read_csv(file2, sep=delimiter)
    elif isinstance(file2, pd.DataFrame):
        df2 = file2
    else:
        raise ValueError("file2 must be a filepath or a DataFrame")
        
    # Rename the merge column in file2 if necessary
    #if file2_merge_col != file1_merge_col:
    #    df2.rename(columns={file2_merge_col: file1_merge_col}, inplace=True)

    # Check if specific columns are to be included from the second file
    if include_columns != "All":
        if not set(include_columns).issubset(df2.columns):
            missing_cols = set(include_columns) - set(df2.columns)
            raise ValueError(f"The following columns are specified in include_columns but do not exist in the second dataset: {missing_cols}")
        columns_to_include = [file1_merge_col] + include_columns
        df2 = df2[columns_to_include]

    # Merge the dataframes on the specified column
    #merged_df = pd.merge(df1, df2, on=file1_merge_col, how=join_type)
    merged_df = pd.merge(df1, df2, left_on=file1_merge_col, right_on=file2_merge_col, how=join_type)
    
    # Drop specified columns if drop_cols is not None or empty
    if drop_cols and len(drop_cols) > 0:
        merged_df = merged_df.drop(drop_cols, axis=1)

    return merged_df

###############################################################################
# Example usage:
# Merging and including all columns from the second file (default behavior)
merged_data_all = merge_data(file1=up_spneumo_proteomes, 
                                  file2=exclusion_reasons,
                                  file1_merge_col='exclusion_id', 
                                  file2_merge_col='id', 
                                  include_columns="All",
                                  join_type='left',
                                  delimiter ='\t',
                                  drop_cols = ['id'])

upid_nunique_source = merged_data_all['upid'].nunique()
print(f"\nLength of source data: {len(merged_data_all)}")
print(f"Number of unique upids in source data: {upid_nunique_source}")

print(merged_data_all.head())


#a = merged_data_all[merged_data_all['gc_set_acc'].isin(['GCA_001255215.2', 'GCA_001255215.1'])]
# Merging and specifying only certain columns to include from the second file
#merged_data_specific = merge_tsv_files(file1_path=up_spneumo_proteomes, 
#                                  file2_path=exclusion_reasons,
#                                  include_columns=['id_description'])

#print("\nMerged with specific columns:")
#print(merged_data_specific.head())

###############################################################################
# checkM: See next!
# FIXME: later!
checkM_filepath = basedir + "/CheckM_report_prokaryotes.txt"
#checkM_ncbi = pd.read_csv(basedir + "/CheckM_report_prokaryotes.txt", sep = "\t")
checkM_ncbi = pd.read_csv(checkM_filepath, sep = "\t")
upcheckM1= merged_data_all[merged_data_all['gc_set_acc'].isin(checkM_ncbi['#genbank-accession'])]
upNOTcheckM1 = merged_data_all[~merged_data_all['gc_set_acc'].isin(checkM_ncbi['#genbank-accession'])]

checkM_ncbi.columns

checkM_ncbi['gc_set_acc']
checkM_ncbi['#genbank-accession']

df_with_checkM = merge_data(file1=merged_data_all,
                                file2=checkM_ncbi,
                                file1_merge_col="gc_set_acc", 
                                file2_merge_col="#genbank-accession", 
                                include_columns=['#genbank-accession', 'refseq-accession', 'checkm-completeness','checkm-contamination'],
                                join_type='left',
                                delimiter = '\t')

df_with_checkM['species_taxid'].value_counts()
df_with_checkM['species-taxid'].value_counts()

df_with_checkM['proteome_taxid'].value_counts()
df_with_checkM['taxid'].value_counts()


# 
# Drop rows where either 'species_taxid' or 'species-taxid' is NaN
small_df = df_with_checkM.dropna(subset=['species_taxid', 'species-taxid', 'proteome_taxid', 'taxid'])

# Perform the comparison and create a new column with the result
small_df['comparison_result'] = small_df['species_taxid'] == small_df['species-taxid']
small_df['comparison_result2'] = small_df['proteome_taxid'] == small_df['taxid']

all_equal = small_df['comparison_result'].all()
print("All non-NaN values are equal:", all_equal)

all_equal2 = small_df['comparison_result2'].all()
print("All non-NaN values are equal:", all_equal2)


print(f"\nLength of checkM_ncbi: {len(checkM_ncbi)}")

print(f"\nLength of available checkM for spneumo acc to GA: {len(upcheckM1)}")
print(f"\nTotal unique UP proteomes with available checkM according to GA: {upcheckM1['upid'].nunique()}") #10161

print(f"\nLength of unavailable checkM for spneumo acc to GA: {len(upNOTcheckM1)}")
print(f"\nTotal unique UP proteomes with UNavailable checkM according to GA: {upNOTcheckM1['upid'].nunique()}") #16588

#checkM_sp = pd.read_csv(basedir + "/CheckM_report_streptococcus_pneumoniae.txt", sep = "\t")
#upcheckM = merged_data_all[merged_data_all['gc_set_acc'].isin(checkM_sp['#genbank-accession'])]
#upNOTcheckM = merged_data_all[~merged_data_all['gc_set_acc'].isin(checkM_sp['#genbank-accession'])]

#(upcheckM1 == upcheckM).all()
#(upNOTcheckM1 == upNOTcheckM).all()

#################
# code for is_effective == t and exclusion_ids 1,7,9,14,29,94,96,99

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
df_unique = df.drop_duplicates(['upid', 'protein_count'])

print(f"\nLength of input data: {len(df)}")
print(f"\nLength of input data unique proteomes: {df['upid'].nunique()}")
print(f"\nCount effective exclusions:\n {df['is_effective'].value_counts()}")
print(f"\nUnique Count effective exclusions:\n {df_unique['is_effective'].value_counts()}")

print(f"\nExcluding proteomes with protein count: {excluded_protein_counts}")
df1 = df[~df['protein_count'].isin(excluded_protein_counts)]
df1_unique = df_unique[~df_unique['protein_count'].isin(excluded_protein_counts)]

print(f"\nLength of unique proteomes in filtered data: {df1['upid'].nunique()}")
print(f"\nLength of df: {len(df1)}")
print(f"\nNo of unique proteomes removed: {df['upid'].nunique() - df1['upid'].nunique()}")
print(f"\nSanity check2 : No of unique proteomes removed: {len(df_unique) - len(df1_unique)}")
a = df_unique['protein_count'].value_counts()
a[0] #76 should match this

print(f"\nExcluding proteomes with assembly level: {assembly_level_to_exclude}")
df2 = df1[~df1['assembly_level'].isin(assembly_level_to_exclude)]
df2_unique =  df1_unique[~df1_unique['assembly_level'].isin(assembly_level_to_exclude)]

print(f"\nLength of df: {len(df2)}")
print(f"\nNo of unique proteomes removed: {df1['upid'].nunique() - df2['upid'].nunique()}")
print(f"\nSanity check2 : No of unique proteomes removed: {len(df1_unique) - len(df2_unique)}")
a2 = df1_unique['assembly_level'].value_counts()
a2['partial']


#4 upids belongig to 10 records have partial
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

# Usage
ids_to_keep = [29, 94, 96, 99]
ids_to_omit = [1, 7, 9, 14, 2]    #26546


# Step 1: Identify upids to remove
grouped = df2.groupby('upid')['exclusion_id'].agg(set)
#grouped2 = df2.groupby('upid')['exclusion_id'].agg(list)
print(f"\nLength of unique proteomes in grouped: {len(grouped)}")

df2.groupby('upid').nunique() # should match

upids_to_remove = grouped[
    grouped.apply(lambda ids: any(i in ids_to_omit for i in ids) and not any(i in ids_to_keep for i in ids))
].index

if len(ids_to_keep) > 0: 
    upids_to_filter = grouped[
        grouped.apply(lambda ids: any(i in ids_to_keep for i in ids) and any(i in ids_to_omit for i in ids))
    ].index

print(f"\nLength of unique proteomes to remove: {(upids_to_remove).nunique()}")
print(f"\nLength of unique proteomes to filter: {upids_to_filter.nunique()}")
#print(f"\nLength of Total unique proteomes to : {upids_to_filter.nunique()}")


upids_to_remove1 = grouped[
    grouped.apply(lambda ids: any(i in ids_to_omit for i in ids))].index
print(f"\nLength of unique proteomes to remove with second criteria: {(upids_to_remove1).nunique()}")


a = upids_to_remove.to_list() + upids_to_filter.to_list()
len(a)
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

print(f"\nNo of unique proteomes in final: {df_grouped['upid'].nunique()}")
print(f"\nNo of unique proteomes removed: {df['upid'].nunique() - df_grouped['upid'].nunique()}")

######

cols_move_end = ['exclusion_id','id_description','is_effective','source']
df_grouped_reordered = move_cols_to_end(df_grouped, cols_move_end)


outfile_df = basedir + "/spneumo_biosample_info_processed.out"
df_grouped_reordered.to_csv(outfile_df, sep="\t", index=False )

outfile_ids = basedir + "/proc-ids.txt"
df_grouped_reordered['upid'].to_csv(outfile_ids, index=False)

#############################

