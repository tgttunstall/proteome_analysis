#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:15:23 2024

@author: tanu
"""
#import time
#import pandas as pd
#import numpy as np
#import os, sys
#from glob import glob
#from tqdm import tqdm
#from Bio import SeqIO # pip install biopython

###############################################################################
#TODO: Add cmd args
#==================
# source dir
#==================
# # ebi-cluster
# # = "/nfs/production/martin/uniprot/users/insana/proteomescomparison"

# # home
# source_dir = "/home/pub/Work/data_arise_proteome"

# # ebi-local machine
# #source_dir = "/home/tunstall/Documents/arise/data_arise_proteome"

# #==================
# # files
# #==================
# # pig proteome file

# ## all .fa files
# proteome_indir = source_dir + "/pig"
# print("\nReading proteome file from:", proteome_indir)

# ## results files: .tsv
# results_indir = source_dir + "/results_pig"
# protein_cluster_infile = results_indir + "/Specie_protein_cluster.tsv" 
# #protein_cluster_infile = results_indir + "/sample_Specie_protein_cluster.tsv"

# #outdir = results_indir
###############################################################################
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Proteome comparison script with configurable input and output paths.")
parser.add_argument("--fasta_dir", type=str, default="/home/pub/Work/data_arise_proteome/pig",
                    help="Directory containing all .fa files (default: /home/pub/Work/data_arise_proteome/pig)")
parser.add_argument("--tsv_file", type=str, default="/home/pub/Work/data_arise_proteome/results_pig/Specie_protein_cluster.tsv",
                    help="Path to the TSV file containing protein clusters (default: Specie_protein_cluster.tsv)")
parser.add_argument("--outfile", type=str, default="/home/pub/Work/data_arise_proteome/results_pig/species_protein_cluster_with_proteome.tsv",
                    help="Path to the output file (default: species_protein_cluster_with_proteome.tsv)")

# Parse arguments
args = parser.parse_args()

# Source directory
source_dir = os.path.dirname(args.fasta_dir)

# Directories and files based on arguments
proteome_indir = args.fasta_dir
protein_cluster_infile = args.tsv_file
output_file = args.outfile

# Output the paths being used
print("\nConfigurations:")
print(f"  Reading proteome files from: {proteome_indir}")
print(f"  Reading protein cluster TSV file from: {protein_cluster_infile}")
print(f"  Writing output to: {output_file}")

