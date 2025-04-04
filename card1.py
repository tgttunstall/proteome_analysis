#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 13:22:38 2025

@author: tanu
"""
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

###############################################################################
# Read data
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/card"

#input_file  = basedir+"/DEL/clusterizes_proteomes2"
#input_file  = basedir+"/Count_Labelled_Species_protein_cluster.tsv"

# proteins
infile_cardIndex = basedir+"/aro_index.tsv"
infile_cardOntology = basedir+"/aro.tsv"
infile_cardUP = basedir+"/CARD-UniProt-Mapping.tsv"

##############################################################################

cardIndex_df = pd.read_csv(infile_cardIndex, sep='\t', header=0 )
cardOntology_df = pd.read_csv(infile_cardOntology, sep='\t', header=0 )
cardUP_df = pd.read_csv(infile_cardUP, sep='\t', header=0 )

cardIndex_df.columns
cardIndex_df.columns = cardIndex_df.columns.str.replace(' ', '_')
cardIndex_df.columns

cardOntology_df.columns
cardOntology_df.rename(columns={'Accession': 'ARO Accession'}, inplace=True)
cardOntology_df.columns = cardOntology_df.columns.str.replace(' ', '_')
cardOntology_df.columns

cardUP_df.columns
cardUP_df.columns = cardUP_df.columns.str.replace(' ', '_')
cardUP_df.columns

##############################################################################
# TEST : Perform the left join
CI_df = pd.merge(cardOntology_df, cardIndex_df, on='ARO_Accession', how='left', indicator='data_source_merge1')
CI_df['data_source_merge1'].value_counts()

duplicates = CI_df[CI_df.duplicated('ARO_Accession', keep=False)]
##############################################################################
def custom_aggregate(series):
    # Convert all items in the series to string and join with a comma
    return ', '.join(series.astype(str))

# Creating an aggregation dictionary dynamically
agg_funcs = {col: custom_aggregate for col in cardIndex_df.columns if col != 'ARO_Accession'}

#------
# MERGE 1: value merging in cardIndex_df
#------
# Group by 'ARO_Accession' and aggregate with custom function
cardIndex_df2 = cardIndex_df.groupby('ARO_Accession').agg(agg_funcs).reset_index()
duplicates1 = cardIndex_df2[cardIndex_df2.duplicated('ARO_Accession', keep=False)]

card_mergeOI_df=pd.merge(cardOntology_df, cardIndex_df2, on='ARO_Accession', how='left', indicator='data_source_merge_ont_ind')
card_mergeOI_df=card_mergeOI_df.drop('data_source_merge_ont_ind', axis=1)
card_mergeOI_df.columns
card_mergeOI_df.rename(columns={'CARD_Short_Name_y': 'CARD_Short_Name_IndexFile',
                                'CARD_Short_Name_x': 'CARD_Short_Name_OntologyFile'}, inplace=True)
card_mergeOI_df.columns

#------
# MERGE 2: OI and CardUP
# want to get everything that is not matched in cardUP 
#------
# Perform an outer join and include an indicator
merged_df = pd.merge(card_mergeOI_df, cardUP_df, on='ARO_Accession', how='outer', indicator='data_source')
merged_df['data_source'].value_counts()

# Filter to get rows that only appear in one of the dataframes 
unmatched_df = merged_df[merged_df['data_source'] != 'both']
unmatched_df['data_source'].value_counts()

# Display the result
print(unmatched_df)

# Rename the values in the custom 'Merge_Indicator' column
unmatched_df['data_source'] = unmatched_df['data_source'].replace({
    'left_only': 'CARD_Ontology_Index',
    'right_only': 'CARD_UP_Mapping',
    'both': 'Both'
})

unmatched_df['data_source'].value_counts()

unmatched_df.columns
unmatched_df['CARD_Short_Name_OntologyFile']==unmatched_df['CARD_Short_Name_IndexFile']
f=unmatched_df['CARD_Short_Name_OntologyFile']==unmatched_df['CARD_Short_Name_IndexFile']
f.value_counts()

# card short names: 51 mismatch
card_short_names_mismatch_df=unmatched_df[unmatched_df['CARD_Short_Name_OntologyFile']!=unmatched_df['CARD_Short_Name_IndexFile']]
outfilename_cname=basedir+"/card_short_names_mismatch_ont_ind.tsv"
card_short_names_mismatch_df.to_csv(outfilename_cname, sep = '\t', index=None)

# card resistance mechanisms: 1994 mismatch
f2=unmatched_df['Resistance_Mechanism_x']==unmatched_df['Resistance_Mechanism_y']
f2.value_counts()

card_res_mismatch_df=unmatched_df[unmatched_df['Resistance_Mechanism_x']!=unmatched_df['Resistance_Mechanism_y']]
outfilename_cres=basedir+"/card_res_mismatch_ont_ind.tsv"
card_res_mismatch_df.to_csv(outfilename_cres, sep = '\t', index=None)

############
cardUP_df.columns
a=unmatched_df.isna().sum()

# drop cols with 1994
#cols_to_delL = ['UPKB', 'CARD_Short_Name', 'RM_ARO_Accession', 'Resistance_Mechanism_y', 'CARD_URL']
# Extract column names where the count of NaNs is length of unmatched_df i.e len(unmatched_df)
cols_to_delL= a[a == len(unmatched_df)].index.tolist()

unmatched_df.drop(cols_to_delL, axis=1, inplace=True)
unmatched_df.columns

unmatched_df.columns = unmatched_df.columns.str.replace('_x', '')
#unmatched_df.columns = unmatched_df.columns.str.replace('_y', '')
unmatched_df.columns

unmatched_df.isna().sum()
len(unmatched_df.columns)

outfilename_unmatched = basedir+"/card_ontology_index_not_mapped_to_UP.tsv"
unmatched_df.to_csv(outfilename_unmatched, sep = '\t', index=None)
print(f"Output File being written to: {outfilename_unmatched}")
print(f"Shape of output file: {unmatched_df.shape}")

