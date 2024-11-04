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
def find_proteome_id(protein_name, proteomeDict):
    """
    Returns the proteome ID from proteomeD for a given protein name,
    or 'XXX' if the protein is not found in proteomeD.
    """
    return proteomeDict.get(protein_name, "XXX")

#     if protein_name in proteomeDict:
#         proteomeDict[protein_name] += f"_{proteome_id}"
#     else:
#         proteomeDict[protein_name] = proteome_id

###############################################################################
# # Attempt 1
# def find_proteome_id(protein_name, proteomeDict):
#     """
#     Function to find corresponding proteome ID for a given protein ID
    
#     Parameters
#     ----------
#     protein_name : string
#         comes from a tsv file of the format:
#             ENSSSCP00000055324|661	ENSSSCP00000055324|661
#             ENSSSCP00000080445|426	ENSSSCP00015043382|461
#     proteinDict : dict
#         contains proteome ids as keys: list of corresponding proteins as values

#     Returns proteome ID or an empty string if not found.
#     -------
#     """
#     for proteome_id, protein_id in proteomeDict.items():
#       if protein_name in protein_id:
#           #print(proteome_id)
#           return proteome_id
#     return "XXX"

#first2pairs = {k: proteomeD[k] for k in list(proteomeD)[:2]}

#TODO: add unit tests here
# test
#find_proteome_id("ENSSSCP00050041795|167", proteomeD)
#find_proteome_id("ENSSSCP00000055324|661", proteomeD)
#find_proteome_id("", proteomeD)

#for proteome_id, protein_id in proteomeD.items():
    #print (proteome_id, (protein_id[0]), print(type(protein_id[0])))
#    if "ENSSSCP00050041795|167" in protein_id:
#        print(proteome_id)

###############################################################################
