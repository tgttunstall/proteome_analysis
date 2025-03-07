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
###############################################################################
def merge_data(file1, file2, 
               file1_merge_col='exclusion_id', 
               file2_merge_col='id', 
               join_type='left',
               delimiter="\t",
               drop_cols=['id']):
    """
    Merges two datasets based on a merge column. Includes all columns from both datasets.
    The type of join performed can be specified to control how rows from each dataset are included in the result.
    Each dataset can either be a DataFrame or a file path to a TSV file.

    Args:
    - file1 (str or pd.DataFrame): First dataset or path to the TSV file.
    - file2 (str or pd.DataFrame): Second dataset or path to the TSV file, containing additional data.
    - file1_merge_col (str): Column name in the first dataset to merge on. Default is 'exclusion_id'.
    - file2_merge_col (str): Column name in the second dataset to merge on. Default is 'id'.
    - join_type (str): Type of join to perform during the merge. Options include 'left', 'right', 'outer', 'inner'. Default is 'left'.
    - delimiter (str): Delimiter used in the file, default is tab ('\t').
    - drop_cols (list): Columns to drop after the merge.

    Returns:
    - pd.DataFrame: A DataFrame with the merged content.

    Raises:
    - ValueError: If there are any issues loading the data.
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
        
    # Merge the dataframes on the specified columns
    merged_df = pd.merge(df1, df2, left_on=file1_merge_col, right_on=file2_merge_col, how=join_type)
    
    # Drop specified columns if drop_cols is not None or empty
    if drop_cols:
        # Check if the columns to be dropped exist in the merged DataFrame
        missing_drop_cols = [col for col in drop_cols if col not in merged_df.columns]
        if missing_drop_cols:
            raise ValueError(f"The following columns to drop do not exist in the DataFrame: {missing_drop_cols}")
        merged_df = merged_df.drop(drop_cols, axis=1)
        
    #TODO:
    # Add sanity check to check the no. of rows and cols
    
    return merged_df

###############################################################################

def process_up_exclusions(df,
                           exclusion_id_colname='exclusion_id',
                           grouping_colname = 'upid',
                           protein_count_colname = 'protein_count',
                           colnames_for_value_merging=[],
                           ids_to_keep=[], 
                           ids_to_omit=[], 
                           excluded_protein_counts=[0], 
                           assembly_level_to_exclude=['partial', 'full'],
                           assembly_colname='assembly_level',
                           reorder_cols=True
                           ):
    
    """
    Processes exclusion reasons from a dataset, handling cases where multiple exclusion reasons
    may exist for the same UPID. The function applies logic to keep or omit records based on user-defined criteria,
    and manages the aggregation of values when multiple records are merged into a single entry.

    Args:
    - df (pd.DataFrame): DataFrame containing the protein data.
    - exclusion_id_colname (str): The name of the column containing exclusion IDs. Default is 'exclusion_id'.
    - grouping_colname (str): The name of the column to group data by, typically UPID. Default is 'upid'.
    - protein_count_colname (str): The column name where protein counts are stored. Default is 'protein_count'.
    - colnames_for_value_merging (list): List of column names whose values should be merged during grouping.
    - ids_to_keep (list): List of exclusion IDs that, if found within a group, will prevent records in that group from being omitted.
    - ids_to_omit (list): List of exclusion IDs that, if found within a group, will cause records in that group to be omitted.
        NB: If a record occurs in both ids_to_keep AND ids_to_omit, the record for the ids_to_keep will take precedence. If there are multiple, only one will be kept.
        This is best left empty to only basically just let the ids_to_omit filter out the bad records! if you still need this complex filtering, it is here to do that:-)
    - excluded_protein_counts (list): Protein counts that should lead to the exclusion of the record.
    - assembly_level_to_exclude (list): List of assembly levels that should lead to the exclusion of the record.
    - reorder_cols (bool): If True, reorders columns after processing, placing 'colnames_for_value_merging' at the end.

    Returns:
    - pd.DataFrame: The processed DataFrame with UPIDs and their exclusion reasons collapsed, filtered by the provided criteria.
    """
    
    ###########################
    # internal functions used
    def move_cols_to_end(df, cols):
        return df[[x for x in df.columns if not x in cols] + cols]
    
    # Define aggregation functions for each column
    
    # This is to check if protein_count during grouping is a single value.
    def check_single_value(series):
        unique_values = set(series)
        if len(unique_values) == 1:
            return next(iter(unique_values))  # Convert the single element set back to an integer
        else:
            raise ValueError("More than one unique value found in protein_count")

    # Update aggregation functions dictionary
    agg_funcs = {
        #col: (lambda x: set(x) if col == 'protein_count' else ', '.join(str(i) for i in x.unique()))
        #col: (lambda x: check_single_value if col == protein_count_colname else ', '.join(str(i) for i in x.unique()))
        col: (lambda x, col=col: check_single_value(x) if col == protein_count_colname else ', '.join(str(i) for i in x.unique()))
        for col in colnames_for_value_merging+[protein_count_colname]
    }

            
    ################################################
    print(f"\nLength of input data: {len(df)}")
    print(f"\nLength of input data unique proteomes: {df['upid'].nunique()}")

    # STEP 0:
    print("\nSTEP 0")
    print(f"|--Excluding proteomes with protein count: {excluded_protein_counts}")
    df1 = df[~df[protein_count_colname].isin(excluded_protein_counts)]
    print(f"|--Length of data after excluding such proteomes: {len(df1)}")
    print(f"|--Length of unique proteomes in data: {df1[grouping_colname].nunique()}")
    print(f"|--No. of unique proteomes removed: {df[grouping_colname].nunique() - df1[grouping_colname].nunique()}")    
    p_zero_count = df[df[protein_count_colname].isin(excluded_protein_counts)][grouping_colname].unique().tolist()

    # STEP 1: 
    print("\nSTEP 1")
    print(f"|--Excluding proteomes with assembly level: {assembly_level_to_exclude}")
    df2 = df1[~df1[assembly_colname].isin(assembly_level_to_exclude)]
    print(f"|--Length of data after excluding such proteomes: {len(df2)}")
    print(f"|--Length of unique proteomes in data: {df2[grouping_colname].nunique()}")
    print(f"|--No. of unique proteomes removed: {df1[grouping_colname].nunique() - df2[grouping_colname].nunique()}")
    p_partial_count = df1[df1[assembly_colname].isin(assembly_level_to_exclude)][grouping_colname].unique().tolist()

    # STEP 2: Grouping by 
    print("\nSTEP 2")
    print(f"|--Grouping data by: {grouping_colname}")
    print(f"|----and aggregating values for: {exclusion_id_colname}")
    #grouped = df2.groupby('upid')['exclusion_id'].agg(set) #or  list 
    grouped = df2.groupby(grouping_colname)[exclusion_id_colname].agg(set)
    print(f"|--Length of unique proteomes in grouped data: {len(grouped)}")

    # STEP 3: Gathering upids to remove (and filter if any)
    print("\nSTEP 3")
    print(f"|--Gathering upids to remove based on {exclusion_id_colname}: {ids_to_omit}")
    
    if len(ids_to_keep) > 0:
        print(f"|--ids_to_keep provided : {ids_to_keep}, hence complex filtering will be applied")
        upids_to_remove = grouped[
            grouped.apply(lambda ids: any(i in ids_to_omit for i in ids) and not any(i in ids_to_keep for i in ids))].index
        print(f"\nNo. of {grouping_colname} to remove: {upids_to_remove.nunique()}")
        
        upids_to_filter = grouped[
            grouped.apply(lambda ids: any(i in ids_to_keep for i in ids) and any(i in ids_to_omit for i in ids))].index
        print(f"\nNo. of {grouping_colname} to filter: {upids_to_filter.nunique()}")

        # Step 3a: Filter DataFrame on ids_to_omit
        df_filtered = df2[~df2[grouping_colname].isin(upids_to_remove)]  # Remove unwanted upids
        
        # Step 3b: Filter DataFrame on ids_to_keep
        df_filtered = df_filtered[
            df_filtered[grouping_colname].isin(upids_to_filter) & df_filtered[exclusion_id_colname].isin(ids_to_keep) | 
            ~df_filtered[grouping_colname].isin(upids_to_filter)
        ]        
        
        print(f"|--Length of data after removing and filtering {grouping_colname}: {(len(df_filtered))}")
        print(f"|--Length of unique proteomes in the data after removing and filtering {grouping_colname}: {df_filtered[grouping_colname].nunique()}")

    else:
        print(f"|--Just filtering data based on {exclusion_id_colname}: {ids_to_omit}")
        upids_to_remove = grouped[
            grouped.apply(lambda ids: any(i in ids_to_omit for i in ids))].index
        print(f"|--No. of {grouping_colname} to remove: {upids_to_remove.nunique()}")

        # Step 3a: Filter DataFrame on ids_to_omit
        df_filtered = df2[~df2[grouping_colname].isin(upids_to_remove)]  # Remove unwanted upids
        print(f"|--Length of data after removing such {grouping_colname}: {(len(df_filtered))}")
        

    f=list(upids_to_remove)
    df_excluded = df[df[grouping_colname].isin(f+p_zero_count+p_partial_count)]  # Remove unwanted upids

    # STEP 4: Applying groupby with custom aggregation functions
    print("\nSTEP 4")
    print(f"|--Applyinn groupby and aggregation of values for columns.\nTotal columns to process: {len(df_filtered.columns)}")
    print(F"\n|--Grouping column: {grouping_colname}")
    print(f"\n|---For {len(colnames_for_value_merging)} columns: Appplying custom aggregation type 'join' values separated by ','. \nThese are {colnames_for_value_merging}")
    print(f"\n|---- For column name {protein_count_colname}: Appplying custom aggregation type 'set' and then checking whether it contains a single value")

    # For all other columns, use 'first' as the aggregation function
    #other_cols = df_filtered.columns.difference(['upid', 'protein_count', 'exclusion_id']).tolist()
    other_cols = df_filtered.columns.difference([grouping_colname]).tolist()

    other_colnames = df_filtered.columns.difference([grouping_colname, protein_count_colname] + colnames_for_value_merging)
    print(f"\n|----For all other {len(other_colnames)} columns: Appplying aggrgation type 'first'. \nThese are {list(other_colnames)}")

    # Default aggregation function for other columns: first
    for col in other_cols:
        #print (col)
        if col not in agg_funcs:  # Ensure we don't overwrite custom columns
            #print(col)
            agg_funcs[col] = 'first'
    
    #df_grouped = df_filtered.groupby(['upid', 'protein_count']).agg(agg_funcs).reset_index()
    df_grouped = df_filtered.groupby([grouping_colname]).agg(agg_funcs).reset_index()
    print(f"\nLength of unique proteomes in grouped and aggregated data: {(df_grouped[grouping_colname]).nunique()}")

    # Sanity check
    print("\nSTEP 5")
    print(f"\Quick sanity check to see if there are any duplicate {grouping_colname}")
    duplicates = df_grouped[df_grouped.duplicated(subset=['upid', 'protein_count'], keep=False)]
    if not duplicates.empty:
        print("\nFAILURE: Duplicates found in grouped data.")
        raise ValueError("Duplicate records found after grouping and aggregation.")
    else:
        print("\nSUCCESS: Data grouped and aggregated successfully.")
        
    if reorder_cols:
        print("\nSTEP 5a")
        print("|--Reordering columns...")
        print(f"|--Placing {colnames_for_value_merging} at the end...")
        df_grouped = move_cols_to_end(df_grouped, colnames_for_value_merging) 

    print("\nQuick QC")
    print(f"|--No of unique proteomes in final data: {df_grouped[grouping_colname].nunique()}")
    print(f"|--No. of unique proteomes removed: {df[grouping_colname].nunique() - df_grouped[grouping_colname].nunique()}")
    print(f"|--Stats for protein length in the proteomes:\n {df_grouped[protein_count_colname].describe()}")

    #return df_grouped
    return df_grouped, df_excluded

###############################################################################




