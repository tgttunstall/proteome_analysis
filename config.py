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
# ebi-cluster
# = "/nfs/production/martin/uniprot/users/insana/proteomescomparison"

# home
source_dir = "/home/pub/Work/data_arise_proteome"

# ebi-local machine
#source_dir = "/home/tunstall/Documents/arise/data_arise_proteome"

#==================
# files
#==================
# pig proteome file

## all .fa files
proteome_indir = source_dir + "/pig"
print("\nReading proteome file from:", proteome_indir)

## results files: .tsv
results_indir = source_dir + "/results_pig"
protein_cluster_infile = results_indir + "/Specie_protein_cluster.tsv" 
#protein_cluster_infile = results_indir + "/sample_Specie_protein_cluster.tsv"

#outdir = results_indir
###############################################################################
