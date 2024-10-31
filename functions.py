#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:20:14 2024

@author: tanu
"""

def find_proteome_id(protein_name, proteomeDict):
    """
    Function to find corresponding proteome ID for a given protein ID
    
    Parameters
    ----------
    protein_name : string
        comes from a tsv file of the format:
            ENSSSCP00000055324|661	ENSSSCP00000055324|661
            ENSSSCP00000080445|426	ENSSSCP00015043382|461
    proteinDict : dict
        contains protome ids as keys: list of corresponding proteins as values

    Returns proteome ID or an empty string if not found.
    -------
    """
    for proteome_id, protein_id in proteomeDict.items():
      if protein_name in protein_id:
          #print(proteome_id)
          return proteome_id
    return "XXX"

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