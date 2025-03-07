#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 16:11:58 2025

@author: tanu
"""
from processing_functions.py import merge_data, process_up_exclusions

# TODO: add cmd args to call these files

# Read data
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
up_spneumo_proteomes = basedir + "/spneumo_biosample_info.out"
exclusion_reasons = basedir + "/exclusion_reasons.tsv"
checkm = pd.read_csv(basedir + "/CheckM_report_prokaryotes.txt", sep="\t")
################################################################################################################################################################################################################################################################################################################################
a = merge_data(file1=up_spneumo_proteomes,
              file2=exclusion_reasons,
              file1_merge_col='exclusion_id', 
              file2_merge_col='id', 
              join_type='left',
              delimiter ='\t',
              drop_cols = ['id'])

a.shape
###############################################################################
checkm.shape
checkm.columns
a2 = merge_data(file1=a,
              file2=checkm,
              file1_merge_col='gc_set_acc', 
              file2_merge_col='#genbank-accession', 
              join_type='left',
              delimiter ='\t',
              drop_cols = ['#genbank-accession', 'taxid', 'species-taxid', 'organism-name', 'species-name', 'assembly-name', 'checkm-marker-set'])
a2.shape
a2[['refseq-accession','checkm-completeness', 'checkm-contamination']].isna().sum()
a2[['refseq-accession','checkm-completeness', 'checkm-contamination']].count()

###############################################################################
#ids_to_keep = [29, 94, 96, 99]
processed_df = process_up_exclusions(df=a2,
                      exclusion_id_colname='exclusion_id',
                      grouping_colname = 'upid',
                      protein_count_colname = 'protein_count',
                      colnames_for_value_merging=['exclusion_id', 'id_description', 'is_effective', 'source'],
                      ids_to_keep=[], 
                      ids_to_omit=[1, 7, 9, 14, 2], 
                      excluded_protein_counts=[0], 
                      assembly_level_to_exclude=['partial'],
                      reorder_cols=True)

# result of call to this function
outfile_df = basedir + "/spneumo_biosample_info_processed.out"
print(f"\nWriting file: {outfile_df}")
processed_df.to_csv(outfile_df, sep="\t", index=False )

outfile_ids = basedir + "/proc-ids.txt"
print(f"\nWriting ID file: {outfile_ids}")
processed_df['upid'].to_csv(outfile_ids, index=False)
###############################################################################