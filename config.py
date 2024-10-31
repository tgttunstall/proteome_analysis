#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:15:23 2024

@author: tanu
"""

import pandas as pd
import numpy as np
import os, sys
from glob import glob
from Bio import SeqIO # pip install biopython


################################
# work cluster
#source_dir = "nfs/production/martin/uniprot/users/insana/proteomescomparison"

# pig proteome file

# all .fa files from 
#proteome_indir = source_dir + "/pig"
#print("\nReading proteome file from:", proteome_indir)

# results files:
#results_indir = source_dir + "/results_pig"
#protein_cluster_infile = results_indir + "/Specie_protein_cluster.tsv"
#protein_cluster_infile = results_indir + "/sample_Specie_protein_cluster.tsv"

#print("\nReading protein cluster file:", protein_cluster_infile)
################################
# home data store
source_dir = "/home/pub/Work/data_arise_proteome"

# pig proteome file

# all .fa files from 
proteome_indir = source_dir + "/pig"
print("\nReading proteome file from:", proteome_indir)

# results files:
results_indir = source_dir + "/results_pig"
#protein_cluster_infile = results_indir + "/Specie_protein_cluster.tsv"
protein_cluster_infile = results_indir + "/sample_Specie_protein_cluster.tsv"

print("\nReading protein cluster file:", protein_cluster_infile)