#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:17:53 2024

@author: tanu
"""
import time
import pandas as pd
import numpy as np
import os, sys
from glob import glob
from tqdm import tqdm
from Bio import SeqIO # pip install biopython
import pprint as pp
from memory_profiler import memory_usage #pip install ,update yml?

# Enable tqdm for pandas
tqdm.pandas()

# local imports
# Add the directory containing config.py and functions.py to the Python path
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# home
sys.path.append('/home/tanu/git/arise_proteome')

# ebi-local
#sys.path.append('/home/tunstall/git/arise_proteome')

# ebi-cluster


###############################################################################
from config import *  # imports config.py as a module
from functions import *  # imports functions.py as a module

###############################################################################
# Reversing to the new structure:
#protein_to_proteome = {}

#for proteome_id, proteins in proteomeD.items():
#    for protein in proteins:
#        protein_to_proteome[protein] = proteome_id
#627957        
#============
# Stage 1: Create proteome dict
#============
# Start the timer
start_time = time.time()

# Directory containing .fa files
print(f"\nReading files from {proteome_indir}")

# Initialize the protein-to-proteome dictionary as proteomeD
proteomeD = {}

# Find all .fa files
fa_files = glob(os.path.join(proteome_indir, "*.fa"))

# Check if files were found
print(f"Number of .fa files found: {len(fa_files)}")
if len(fa_files) == 0:
    print("No .fa files found. Please check the directory path.")

# Populate proteomeD directly from the .fa files with progress tracking
for file in tqdm(fa_files, desc="Processing .fa files"):
    # Extract proteome ID from the filename (e.g., "proteome_35497" from "proteome_35497.fa")
    proteome_id = os.path.splitext(os.path.basename(file))[0]

    # Read each protein entry in the .fa file and add it to proteomeD
    for record in SeqIO.parse(file, "fasta"):
        #print(f"Processing protein: {record.id} in {proteome_id}")  # Debugging output
        
        #proteomeD[record.id] = proteome_id
        # if record.id in proteomeD:
        #     proteomeD[record.id] += f"_{proteome_id}"
        # else:
        #     proteomeD[record.id] = proteome_id

        if record.id not in proteomeD:
            proteomeD[record.id] = proteome_id        
        else:     
            proteomeD[record.id] += f"_{proteome_id}"

print(proteomeD['ENSSSCP00015020973|188'])

# End the timer and calculate elapsed time
end_time = time.time()
elapsed_time_s1 = end_time - start_time

print(f"Dictionary creation completed in {elapsed_time_s1 // 3600:.0f} hours, "
      f"{(elapsed_time_s1 % 3600) // 60:.0f} minutes, and {elapsed_time_s1 % 60:.2f} seconds.")

# Print sample dictionary content to confirm entries
#print("Sample dictionary entries:", list(proteomeD.items())[627950:627957])
#print("Sample dictionary entries:", list(proteomeD.items())[:5])

print("Sample dictionary entries:")
pp.pprint(list(proteomeD.items())[:5])

#============
# Stage 2: assign proteome id to protein names
#find_proteome_id()
#============

# Read in the TSV file with protein pairs
pclustersDF_original = pd.read_csv(protein_cluster_infile, sep='\t', usecols=[0, 1], names=['protein1_id', 'protein2_id'])

pclustersDF = pclustersDF_original.drop_duplicates()
print(f"\nLength of DF of unique entries:{len(pclustersDF)}")

# Start timing
start_time = time.time()

# Step 1: Add a numeric counter as a new 'cluster_id' column
pclustersDF['cluster_id'] = pclustersDF.groupby('protein1_id').ngroup()

# Step 2: Add a representative column to mark the first instance of each cluster
# We mark the first appearance of each 'protein1_id' as the representative with an asterisk
pclustersDF['representative'] = ''
pclustersDF.loc[pclustersDF.duplicated('protein1_id', keep='first') == False, 'representative'] = '*'

# Step 3: Add 'protein2_id_proteome' column
pclustersDF['protein2_id_proteome'] = pclustersDF['protein2_id'].progress_apply(
    #lambda name: f"{find_proteome_id(name, proteomeD)}_{name}")
    lambda name: f"{find_proteome_id(name, proteomeD)}")

# # without function:
# pclustersDF['protein2_id_proteome'] = pclustersDF['protein2_id'].progress_apply(
#     lambda name: f"{proteomeD.get(name, 'XXX')}_{name}")
#    lambda name: f"{find_proteome_id(name, proteomeD)}")

# Sort the DataFrame by 'cluster_id' in ascending order without resetting the index
pclustersDF = pclustersDF.sort_values(by='cluster_id')

# End timing
end_time = time.time()
elapsed_time_s2 = end_time - start_time

# Output results
print(f"Mapping completed in {elapsed_time_s2 // 3600:.0f} hours, "
      f"{(elapsed_time_s2 % 3600) // 60:.0f} minutes, and {elapsed_time_s2 % 60:.2f} seconds.")
print("Sample mapped entries:", pclustersDF[['cluster_id', 'protein1_id', 'protein2_id', 'representative', 'protein2_id_proteome']].head())


#a  = pclustersDF.iloc[[680512,680513]]

b = pclustersDF[(pclustersDF['cluster_id'] == 301) & (pclustersDF['protein2_id_proteome'].str.contains("ENSSSCP00040031975"))]

pclustersDF.shape

fd_diff = pclustersDF.shape[0] - len(proteomeD)
print(f"\nDiff between tsv file and dict length: {fd_diff} \
      \nLength of TSV file: {pclustersDF.shape[0]} \
      \nLength of proteomeD: {len(proteomeD)}")


duplicateRowsDF = pclustersDF[pclustersDF.duplicated()]
print(f"\nDuplicated rows in pclustersDF: {duplicateRowsDF.shape[0]}")

predicted_nunique_entries = pclustersDF.shape[0] - duplicateRowsDF.shape[0]
print(f"\nPredicted no. of unique entries:{predicted_nunique_entries}")

pclustersDF_unique = pclustersDF.drop_duplicates()
print(f"\nLength of DF of unique entries:{len(pclustersDF_unique)}")

#TODO
# add sanity check here


##############################################################################
#=========================
# writing the updated file
#=========================
# Save the updated DataFrame with proteome ID added
output_clusterDF_file = results_indir + '/species_protein_cluster_with_proteomes_v3.tsv'
cols_to_output = ['cluster_id', 
                  #'protein1_id', 
                  'protein2_id', 
                  'representative', 
                  'protein2_id_proteome']

print(f"writing file: {output_clusterDF_file}\nDim of df: {pclustersDF[cols_to_output].shape}")
pclustersDF[cols_to_output].to_csv(output_clusterDF_file, sep='\t', index=True)
print(f"Updated file saved as {output_clusterDF_file}")

###############################################################################
# FIXME: assumed, that each protein belongs to a unique proteome!
# ISSUE: one to many mapping b/w protein and proteome

#grep -r ENSSSCP00040031975 *.fa
#proteome_4698268.fa:>ENSSSCP00040031975|188
#proteome_4698920.fa:>ENSSSCP00040031975|188
#proteome_4698921.fa:>ENSSSCP00040031975|188
#proteome_4698924.fa:>ENSSSCP00040031975|188

#grep -r ENSSSCP00015020973 *.fa
#proteome_4698268.fa:>ENSSSCP00015020973|188
#proteome_4698920.fa:>ENSSSCP00015020973|188
#proteome_4698921.fa:>ENSSSCP00015020973|188
#proteome_4698924.fa:>ENSSSCP00015020973|188
