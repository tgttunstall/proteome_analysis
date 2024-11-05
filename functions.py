#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:20:14 2024

@author: tanu
"""
###############################################################################
import pandas as pd
import time
from tqdm import tqdm

# Sample proteomeD dictionary for reference (already populated)
# proteomeD = {"ENSSSCP00000055324|661": "proteome_35497", ...}

# Sample data - assuming you've loaded your TSV data into a DataFrame `pclusters`
# pclustersDF = pd.read_csv("protein_cluster.tsv", sep="\t")

# Simplified function to look up the proteome ID
# def find_proteome_id(protein_name, proteomeDict):
#     """
#     Returns the proteome ID from proteomeD for a given protein name,
#     or 'XXX' if the protein is not found in proteomeD.
#     """
#     return proteomeDict.get(protein_name, "XXX")

# #     if protein_name in proteomeDict:
# #         proteomeDict[protein_name] += f"_{proteome_id}"
# #     else:
# #         proteomeDict[protein_name] = proteome_id


def find_proteome_id(protein_name, proteomeDict):
    """
    Returns a string with the format 'protein_name_proteomeID' for a given protein name.
    If the protein is not found in proteomeDict, returns 'protein_name_XXX'.
    """
    proteome_id = proteomeDict.get(protein_name, "XXX")
    return f"{protein_name}_{proteome_id}"


# Define a function to read the file, so memory usage can be tracked
def read_large_tsv(file_path):
    start_time = time.time()
    df = pd.read_csv(file_path, sep='\t')
    end_time = time.time()
    
    elapsed_time = end_time - start_time
    print(f"Time taken to read file: {elapsed_time // 3600:.0f} hours, "
          f"{(elapsed_time % 3600) // 60:.0f} minutes, and {elapsed_time % 60:.2f} seconds.")
    return df
