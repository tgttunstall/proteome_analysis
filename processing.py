#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:17:53 2024

@author: tanu
"""
import pandas as pd
import numpy as np
import os, sys
import glob
from Bio import SeqIO # pip install biopython
#import config.py
################################

# Initialise the dictionary:
pigD = {}

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
    print(f"Reading file: {file}")
    
    # extract proteome_id without the file extension
    proteome_id = os.path.splitext(os.path.basename(file))[0]
    
    # get the headers from the fasta file and store in a list
    headers = [record.id for record in SeqIO.parse(file,"fasta")]
    
    # add the headers list to the proteome_id key in the dict
    pigD[proteome_id] = headers
    
# Print proteome ID and the length of the list of headers for each key
for proteome_id, headers in pigD.items():
    print(f"{proteome_id}: {len(headers)}")