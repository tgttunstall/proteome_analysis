#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import seaborn as sns
#import itertools
#import random

###############################################################################
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
#========
# UP: up db data, checkm, checkm2
#========
# up db data
#up_file = basedir + "/spneumo_biosample_info_v1.out"
#up_file = basedir + "/spneumo_biosample_info_processed.out"
up_file = basedir + "/spneumo_biosample_info.out"

up_all_df = pd.read_csv(up_file, sep="\t")

input_file_up_pid= basedir+"/spneumo_proteomeids.out" 
up_pid_df = pd.read_csv(input_file_up_pid, sep='\t', header=None) #26749,4
up_pid_df.columns=['pid', 'upid', 'gca_set_acc', 'species_name']

# up checkm data
up_checkm = basedir + "/checkm_results_up.tsv"
up_checkm_df = pd.read_csv(up_checkm, sep = "\t")

# up checkm2 data
up_checkm2 = basedir + "/checkm2_results_up.tsv"
up_checkm2_df = pd.read_csv(up_checkm2, sep = "\t")

#========
# ATB
#========
# atb source data
atb_file = basedir + "/atb_pcounts.tsv"
atb_all_df = pd.read_csv(atb_file, sep="\t")

# atb checkm data
atb_checkm = basedir + "/checkm_results_atb.tsv"
atb_checkm_df = pd.read_csv(atb_checkm, sep = "\t")

# atb checkm2 data
atb_checkm2 = basedir + "/checkm2_results_atb.tsv"
atb_checkm2_df = pd.read_csv(atb_checkm2, sep = "\t")

# busco data
#atb_busco = basedir + "/busco_results_atb.tsv"
#atb_busco_df = pd.read_csv(atb_busco, sep = "\t")

###############################################################################
# Sanity check: Filter protein_count>0 and drop duplicates
#========
# UP Data Prep
#========
# First we need to merge up_all_df with pid
pid_df=up_pid_df[['pid', 'upid']]
pid_df['Bin_Id'] = 'proteome_' + pid_df['pid'].astype(str)
pid_df2=pid_df[['Bin_Id', 'upid']]

up_all_pid_df = pd.merge(up_all_df, pid_df2, on='upid', how='left')
# add proteome_pid
up_all_pid_df.columns

up_df = up_all_pid_df[up_all_pid_df['protein_count']> 0]
# Drop duplicates
up_dfU_biosample = up_df.drop_duplicates(subset=['biosample'], keep = 'first')
up_dfU = up_df.drop_duplicates(subset=['upid'], keep = 'first') #26673, 21  [vs 26674 in checkm and checkm2]
up_dfU.columns

# FIXME: 26673 vs 26674 
uL=list(up_dfU['Bin_Id'])
len(uL)

ucmL=list(up_checkm_df['Bin_Id'])
len(ucmL)
a = [i for i, j in zip(set(uL), set(ucmL)) if i != j]
print(a) # proteome_4225276

up_dfU.shape # 26673,21

# Merge dfs for plotting
# checkM and checkM2
up_checkm_df.columns
cm_cols_to_extract=['Bin_Id', 
                 #'Marker_lineage', 
                 #'genomes', 
                 #'markers', 
                 #'marker_sets', '0', '1', '2', '3', '4', '5+', 
                 'Completeness', 
                 'Contamination',
                 'Strain_heterogeneity'
                 #, 'Unnamed: 14'
                 ]
up_cm = up_checkm_df[cm_cols_to_extract]
up_cm.rename({'Completeness' : 'Completeness_CM', 
                'Contamination': 'Contamination_CM',
                'Strain_heterogeneity': 'Strain_heterogeneity_CM'
                }, axis = 1, inplace=True)

up_cm.columns
print(len(up_cm.columns))

up_checkm2_df.columns
cm2_cols_to_extract=['Name', 
                     'Completeness', 
                     'Contamination', 
                     'Completeness_Model_Used'
                     #, 'Additional_Notes'
                     ]

up_cm2 = up_checkm2_df[cm2_cols_to_extract]
up_cm2.rename({'Name': 'Bin_Id',
                'Completeness' : 'Completeness_CM2', 
                'Contamination': 'Contamination_CM2',
                }, axis = 1, inplace=True)
up_cm2.columns
print(len(up_cm2.columns))

up_df_scores = pd.merge(up_cm, up_cm2, on = 'Bin_Id' , how='outer')
print(up_df_scores.head())
print(up_df_scores.columns)
print(up_dfU.columns)

subset_cols_up = ['upid', 
                 # 'proteome_taxid', 
                 # 'species_taxid', 
                  'gc_set_acc', 
                  'biosample',
                  'project_acc', 
                  'protein_count', 
                  'assembly_level',
                  'complete_combined_score',
                  'complete_single_score',
                  'complete_duplicate_score', 
                  'fragment_score', 
                  'missing_score',
                  #'completeness_description', 
                  #'is_representative', 'is_reference', 'is_excluded', 'exclusion_id', 'exclusion_immunity', 'is_redundant',
                  'Bin_Id']

up_dfU_subset = up_dfU[subset_cols_up] #26673,12

up_df_plot = pd.merge(up_dfU_subset, up_df_scores, on='Bin_Id', how = 'outer')
print(up_df_plot.columns)
print(up_df_plot.shape) #26674,18

#===============
# ATB Data Prep
#===============
atb_df = atb_all_df[atb_all_df['protein_count'] > 0]
# Drop duplicates
atb_dfU = atb_df.drop_duplicates(subset=['biosample'], keep = 'first') # 119701,2

# Merge dfs for plotting
# checkM and checkM2
atb_checkm_df.columns
cm_cols_to_extract=['Bin_Id', 
                 #'Marker_lineage', 
                 #'genomes', 
                 #'markers', 
                 #'marker_sets', '0', '1', '2', '3', '4', '5+', 
                 'Completeness', 
                 'Contamination',
                 'Strain_heterogeneity'
                 #, 'Unnamed: 14'
                 ]
atb_cm = atb_checkm_df[cm_cols_to_extract]
atb_cm.rename({'Completeness' : 'Completeness_CM', 
                'Contamination': 'Contamination_CM',
                'Strain_heterogeneity': 'Strain_heterogeneity_CM'
                }, axis = 1, inplace=True)

atb_cm.columns
print(len(atb_cm.columns))

atb_checkm2_df.columns
cm2_cols_to_extract=['Name', 
                     'Completeness', 
                     'Contamination', 
                     'Completeness_Model_Used'
                     #, 'Additional_Notes'
                     ]

atb_cm2 = atb_checkm2_df[cm2_cols_to_extract]
atb_cm2.rename({'Name': 'Bin_Id',
                'Completeness' : 'Completeness_CM2', 
                'Contamination': 'Contamination_CM2',
                }, axis = 1, inplace=True)
atb_cm2.columns
print(len(atb_cm2.columns))

# Merge 1: ATB CM and CM2
#atb_df_scores = pd.merge(atb_cm, atb_cm2, on = 'Bin_Id' , how='outer')
atb_df_scores_cm = pd.merge(atb_cm, atb_cm2, on = 'Bin_Id' , how='outer')

print(atb_df_scores.head())
print(atb_df_scores.columns)
print(atb_dfU.columns)

# Merge 2: merged cm and cm2 + BUSCO
atb_busco_df.columns
atb_busco_df.shape
print(len(atb_busco_df))
atb_df_scores = pd.merge(atb_busco_df, atb_df_scores_cm, left_on='biosample', right_on='Bin_Id', how='outer')

# Merge 3: atb_dfU and atb_df_scores
atb_df_plot = pd.merge(atb_dfU, atb_df_scores, left_on='biosample', right_on='Bin_Id', how = 'outer')
print(atb_df_plot.columns)

###############################################################################
