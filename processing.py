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
import glob
from Bio import SeqIO # pip install biopython
#import config.py
###############################################################################
#============
# Stage 1: Create proteome dict
#============
# Start the timer
start_time_stage1 = time.time()

# Initialise the dictionary:
proteomeD = {}

# Find all .fa files in the directory
fa_files = glob(os.path.join(proteome_indir, "*.fa"))

# Loop through each file and read the sequences
#for file in fa_files:
#    print(f"Reading file: {file}")
    #for record in SeqIO.parse(file, "fasta"):
    #    print(f"ID: {record.id}")
    #    print(f"Sequence: {record.seq[:50]}...")  # Display first 50 bases only    
    
# process each fasta file where you extract the proteome ID from the filename which is the key
# and all corresponding fasta headers as list values

for file in fa_files:
    #print(f"Reading file: {file}")
    
    # extract proteome_id without the file extension
    proteome_id = os.path.splitext(os.path.basename(file))[0]
    
    # get the headers from the fasta file and store in a list
    pfasta_headers = [record.id for record in SeqIO.parse(file,"fasta")]
    
    # add the headers list to the proteome_id key in the dict
    proteomeD[proteome_id] = pfasta_headers
    
# Print proteome ID and the length of the list of headers for each key
#for proteome_id, pfasta_headers in proteomeD.items():
#    print(f"{proteome_id}: {len(pfasta_headers)}")

# End the timer and calculate elapsed time
end_time_stage1 = time.time()
elapsed_time_stage1 = end_time_stage1 - start_time_stage1

print(f"Mapping completed in {elapsed_time_stage1:.2f} seconds.")

# Convert to hours, minutes, and seconds
hours_s1 = int(elapsed_time_stage1 // 3600)
minutes_s1 = int((elapsed_time_stage1 % 3600) // 60)
seconds_s1 = elapsed_time_stage1 % 60
print(f"Mapping completed in {hours_s1} hours, {minutes_s1} minutes, and {seconds_s1:.2f} seconds.")
###############################################################################    
#============
# Stage 2: assign proteom id to protein names
#find_proteome_id()
#============
# Start the timer
start_time_stage2 = time.time()

print("\nReading protein cluster file:", protein_cluster_infile)
#protein_cluster_file = pd.read_csv(protein_cluster_infile, sep='\t', header=None)
# tscv file has no headers, so giving it names
pclustersDF = pd.read_csv(protein_cluster_infile
                          , sep ='\t'
                          , usecols = [0,1], names=['protein1_id', 'protein2_id'])

#pclustersDF['protein1_id'].value_counts() # may be some optimisation can be done here!
#pclustersDF['protein2_id'].value_counts()
print("\n No. of unique protein_ids in the two columns:\n", pclustersDF.nunique())


# Apply the function find_proteome_id to the 2 columns of the tsv file:
## adding proteome id to column: 'protein1_id
#pclusters['protein1_id_proteome'] = pclusters['protein1_id'].apply(lamda name: f"{find_proteome_id(name, proteomeD)}_{name}" if find_proteome_id(name, proteomeD) else name)

pclustersDF['protein1_id_proteome'] = pclustersDF['protein1_id'].apply(
    lambda name: f"{find_proteome_id(name, proteomeD)}_{name}" 
    if find_proteome_id(name, proteomeD) 
    else name
)

## adding proteome id to column: 'protein2_id'
pclustersDF['protein2_id_proteome'] = pclustersDF['protein2_id'].apply(
    lambda name: f"{find_proteome_id(name, proteomeD)}_{name}" 
    if find_proteome_id(name, proteomeD) 
    else name
)

# End the timer and calculate elapsed time
end_time_stage2 = time.time()
elapsed_time_stage2 = end_time_stage2 - start_time_stage2

print(f"Mapping completed in {elapsed_time_stage2:.2f} seconds.")

# Convert to hours, minutes, and seconds
hours_s2 = int(elapsed_time_stage2 // 3600)
minutes_s2 = int((elapsed_time_stage2 % 3600) // 60)
seconds_s2 = elapsed_time_stage2 % 60
print(f"Mapping completed in {hours_s2} hours, {minutes_s2} minutes, and {seconds_s2:.2f} seconds.")

###############################################################################
#=========================
# writing the updated file
#=========================
# Save the updated DataFrame with proteome ID added
updated_clusterDF = results_indir + '/species_protein_cluster_with_proteomes.tsv'
print(f"writing file: {updated_clusterDF}\nDim of df: {pclustersDF.shape}")

pclustersDF.to_csv(updated_clusterDF, sep='\t', index=False)
print(f"Updated file saved as {updated_clusterDF}")
