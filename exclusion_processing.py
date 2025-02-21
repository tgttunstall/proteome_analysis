#!/bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import seaborn as sns
#import random
import argparse
from tqdm import tqdm

######
# EXCLUSION processing
#if a up proteome is excluded (is_excluded == t), then look for exclusion_immunity. :
#    If exclusion_immunity==t, keep it.
#    If exclusion_immunity==f, move on to exclusion_id.
#    If exclusion_id = 2, then omit # partial
#    If exclusion_id = 1, then omit # genome length too large
#    If exclusion_id = 7, then omit # genome length too small
######

# processing exclusions
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
up_file = basedir + "/spneumo_biosample_info.out"

# Read the data
df_up_original = pd.read_csv(up_file, sep="\t")
#print(f"\nDF shape: {df_up_original.shape}")
print(f"\nLength of UP DF: {len(df_up_original)}")

# Step 0: omit the ones with protein_count == 0
df_up_all = df_up_original[df_up_original['protein_count']> 0]
print(f"\nNo. of of entries with 0 proteins: {len(df_up_original[df_up_original['protein_count']== 0])}")
print(f"\nLength of UP DF after removing proteomes with 0 proteins: {len(df_up_all)}")

#####################
# Filter for non-excluded proteomes directly
non_excluded_proteomes = df_up_all[df_up_all['is_excluded'] == 'f']

# Filter for excluded proteomes
excluded_proteomes = df_up_all[df_up_all['is_excluded'] == 't']

# From excluded, select those with exclusion immunity
immune_proteomes = excluded_proteomes[excluded_proteomes['exclusion_immunity'] == 't']

# From excluded, select those to process further (no immunity)
non_immune_proteomes = excluded_proteomes[excluded_proteomes['exclusion_immunity'] == 'f']

# IDs to omit based on specific exclusion_id values
ids_to_omit = [1, 2, 7]

# Filter non-immune proteomes that should not be omitted
excluded_keep_proteomes = non_immune_proteomes[~non_immune_proteomes['exclusion_id'].isin(ids_to_omit)]

# Combine all valid proteomes: non-excluded, immune, and processed non-immune
df_up_keep = pd.concat([non_excluded_proteomes, immune_proteomes, excluded_keep_proteomes])

# Information on the counts and filtering
print("Original Count:", df_up_all.shape[0])
print("Final Count:", df_up_keep.shape[0])
print("Included from Excluded (Immune and Processed):", len(immune_proteomes) + len(excluded_keep_proteomes))

# Example to check data
print(df_up_keep.head())

## sanity check
df_up_all.groupby(['is_excluded']).size()
df_up_all.groupby(['exclusion_immunity', 'is_excluded']).size()
df_up_all.groupby(['is_excluded','exclusion_id']).size()

# Total original count
original_count = len(df_up_all)

# Count non-excluded directly
non_excluded_count = len(df_up_all[df_up_all['is_excluded'] == 'f'])

# Count immune excluded
immune_excluded_count = len(df_up_all[(df_up_all['is_excluded'] == 't') & (df_up_all['exclusion_immunity'] == 't')])

# Count non-immune excluded that should not be omitted
processed_non_immune_count = len(df_up_all[(df_up_all['is_excluded'] == 't') &
                                           (df_up_all['exclusion_immunity'] == 'f') &
                                           (~df_up_all['exclusion_id'].isin(ids_to_omit))])

# Expected valid proteomes count
expected_final_count = non_excluded_count + immune_excluded_count + processed_non_immune_count

# Final dataset count (assuming it's been processed as described previously)
final_count = len(df_up_keep)  # final_dataset is the resulting DataFrame after filtering

# Print counts for verification
print("Original Count:", original_count)
print("Non-Excluded Count:", non_excluded_count)
print("Immune Excluded Count:", immune_excluded_count)
print("Processed Non-Immune Count:", processed_non_immune_count)
print("Expected Final Count:", expected_final_count)
print("Actual Final Count:", final_count)

# Sanity Check
assert expected_final_count == final_count, "Mismatch in counts, check exclusion logic!"

if expected_final_count == final_count:
    print("Sanity check passed: All counts are as expected!")
else:
    print("Sanity check failed: Counts do not match expected values.")
    
    
    
# check
#	upid
#26283	UP000242121
    
##########
# process all duplicates
#upid_counts = df_up_keep['upid'].value_counts()
#df_up_keep['upid'].nunique()
#a = df_up_keep[df_up_keep['upid'].duplicated() == True]

#upid_counts[upid_counts == 4]
#upid_counts[upid_counts == 3]

#c4 = ['UP000092123', 'UP000046746', 'UP000048098', 'UP000044220','UP000038225' ]
#c4 = ['UP000044220','UP000038225' , 'UP000243565']

#a2 = df_up_keep[df_up_keep['upid'].isin(c4)]

###
# identify duplicates
# Select all columns except 'exclusion_id' for duplicate checking
cols_to_check = df_up_keep.columns.difference(['exclusion_id'])

# Identify duplicates based on 'upid' and all other columns except 'exclusion_id'
df_up_keep['is_duplicate'] = df_up_keep.duplicated(subset=['upid'] + cols_to_check.tolist(), keep=False)

# Extract only the duplicates for further examination or processing
duplicates_df = df_up_keep[df_up_keep['is_duplicate']]

# Optionally, view some of these duplicates to check
print(duplicates_df.head())

# Check the number of identified duplicates
print("Number of identified duplicates:", duplicates_df.shape[0])

##
# Step 1: Count initial total entries
initial_count = len(df_up_keep)

# Step 2: Identify duplicates and count them
# Mark duplicates for the specified columns excluding 'exclusion_id'
cols_to_check = df_up_keep.columns.difference(['exclusion_id']).tolist()
df_up_keep['is_duplicate'] = df_up_keep.duplicated(subset=cols_to_check, keep=False)

# Count all entries marked as duplicates
duplicate_counts = df_up_keep['is_duplicate'].sum()

# Step 3: Process duplicates and aggregate 'exclusion_id' (as described in previous steps)
# Aggregate duplicates
duplicates_df = df_up_keep[df_up_keep['is_duplicate']]
aggregated_duplicates = duplicates_df.groupby(cols_to_check).agg({
    'exclusion_id': lambda ids: ', '.join(ids.astype(str))
}).reset_index()

# Remove original duplicates
df_up_keep_cleaned = df_up_keep[~df_up_keep['is_duplicate']]

# Merge aggregated duplicates back
df_up_keep_final = pd.concat([df_up_keep_cleaned, aggregated_duplicates], ignore_index=True)

# Drop the 'is_duplicate' column as it's no longer needed
df_up_keep_final.drop(columns='is_duplicate', inplace=True)

# Step 4: Count entries after processing
final_count = len(df_up_keep_final)

# Expected count after processing
# We assume each group of duplicates is replaced by one entry
unique_duplicates = len(aggregated_duplicates)
expected_final_count = (initial_count - duplicate_counts) + unique_duplicates

# Print results to verify
print(f"Initial count: {initial_count}")
print(f"Duplicate entries identified: {duplicate_counts}")
print(f"Unique duplicate groups: {unique_duplicates}")
print(f"Final count after processing: {final_count}")
print(f"Expected final count: {expected_final_count}")

# Perform the sanity check
assert final_count == expected_final_count, "Sanity check failed: Final counts do not match expected values."
print("Sanity check passed: All counts are as expected!")

###


a  = df_up_keep_final[df_up_keep_final['upid'].duplicated() == True]

df_up_keep_final['protein_count'].describe()

