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

# local imports
# Add the directory containing config.py and functions.py to the Python path
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))

#sys.path.append('/home/tanu/git/arise_proteome')
sys.path.append('/home/tunstall/git/arise_proteome')

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
        print(f"Processing protein: {record.id} in {proteome_id}")  # Debugging output
        proteomeD[record.id] = proteome_id

# End the timer and calculate elapsed time
end_time = time.time()
elapsed_time_s1 = end_time - start_time

print(f"Dictionary creation completed in {elapsed_time_s1 // 3600:.0f} hours, "
      f"{(elapsed_time_s1 % 3600) // 60:.0f} minutes, and {elapsed_time_s1 % 60:.2f} seconds.")

# Print sample dictionary content to confirm entries
print("Sample dictionary entries:", list(proteomeD.items())[627950:627957])

#============
# Stage 2: assign proteome id to protein names
#find_proteome_id()
#============
# Start the timer
start_time = time.time()

print("\nReading protein cluster file:", protein_cluster_infile)
#protein_cluster_file = pd.read_csv(protein_cluster_infile, sep='\t', header=None)

# tsv file has no headers, so giving it names
pclustersDF = pd.read_csv(protein_cluster_infile
                          , sep ='\t'
                          , usecols = [0,1], names=['protein1_id', 'protein2_id'])

print(pclustersDF.shape)
#TODO: room for optimisation
#pclustersDF['protein1_id'].value_counts() # may be some optimisation can be done here!
#pclustersDF['protein2_id'].value_counts()
print("\n No. of unique protein_ids in the two columns:\n", pclustersDF.nunique())
check1 = pclustersDF.groupby('protein1_id').count()


# Start timing
#start_time = time.time()

# Apply the function to `protein2_id` column with progress tracking
tqdm.pandas(desc="Mapping protein2_id to proteome ID")

## Adding proteome id to columns: 

# FIXME: Adding to column: 'protein1_id'
pclustersDF['protein1_id_proteome'] = pclustersDF['protein1_id'].progress_apply(
    lambda name: f"{find_proteome_id(name, proteomeD)}_{name}"
)

# Adding to column: 'protein2_id'
pclustersDF['protein2_id_proteome'] = pclustersDF['protein2_id'].progress_apply(
    lambda name: f"{find_proteome_id(name, proteomeD)}_{name}"
)

# End timing
end_time = time.time()
elapsed_time_s2 = end_time - start_time

# Output results
print(f"Mapping completed in {elapsed_time_s2 // 3600:.0f} hours, "
      f"{(elapsed_time_s2 % 3600) // 60:.0f} minutes, and {elapsed_time_s2 % 60:.2f} seconds.")
print("Sample mapped entries:", pclustersDF[['protein2_id', 'protein2_id_proteome']].head())

###############################################################################
#=========================
# writing the updated file
#=========================
# Save the updated DataFrame with proteome ID added
updated_clusterDF = results_indir + '/species_protein_cluster_with_proteomes.tsv'
print(f"writing file: {updated_clusterDF}\nDim of df: {pclustersDF.shape}")

pclustersDF.to_csv(updated_clusterDF, sep='\t', index=False)
print(f"Updated file saved as {updated_clusterDF}")

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
