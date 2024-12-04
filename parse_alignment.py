#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:44:27 2024

@author: tunstall
"""

from Bio import Phylo
from Bio import AlignIO
import pandas as pd
import os, sys
import numpy as np

homedir = os.path.expanduser("~")

#######################
#that prints the columns for variations
#
#id               221 223
#35497:ENSCP....  K M
#123:ENSCP...     K N
#2323:ENSCP..     L P
#####################

ALIGN_FORMAT = "fasta"

#aln_fname = "/home/tunstall/git/arise_proteome/11576.aln"
#aln_fname = "/home/tunstall/git/arise_proteome/829.aln"
aln_fname = "/home/pub/Work/data_arise_proteome/testC/results_testC_v3/20613.aln"

align = AlignIO.read(aln_fname, ALIGN_FORMAT)
align[:,0]

###############################################################################
# TODO 1: time this with diff formats: clustal, fasta, 

#TODO 2: you can repeat a timing test using timeit e.g. 1000 times (or 100, or 10 if it's too slow)
# compare the performance of
# align2vardf_for
# vs
# align2vardf_npfor small medium and huge alignment files
# with and without gaps:
    
# tip: use timeit
# https://docs.python.org/3/library/timeit.html
###############################################################################

def align2vardf_for(alignment, nogaps=True):
    #create variation df
    variation_dict = {}

    aln_length = alignment.get_alignment_length()
    
    if nogaps:
        for pos in range(0, aln_length):
            column = alignment[:, pos]

            prev_res = column[0]
            
            if prev_res == "-":
                continue

            found_var = False
            for res in column[1:]:
                if res == "-":
                    found_var = False
                    break
                if res != prev_res:
                    found_var = True
            if found_var:
                print(f"{pos} {column}")
                variation_dict[pos] = list(column)
    else:
        for pos in range(0, aln_length):
            column = alignment[:, pos]

            prev_res = column[0]
            for res in column[1:]:
                if res != prev_res:
                    #print(f"{pos} {column}") #debug
                    variation_dict[pos] = list(column)
                    break

    variation_df = pd.DataFrame.from_dict(variation_dict, orient = "columns")
    variation_df.index = [alignment[i].id for i in range(len(align))]
    return variation_df


def align2vardf_np(alignment, nogaps=True):
    """
    Create a DataFrame containing only the columns from an alignment with variation
    (i.e. where sequence identity is not 100% over all the sequences of the alignment).
    
    Parameters:
    alignment -- Bio.Align object
    nogaps -- bool, if True, excludes columns with gaps
    
    Returns:
    A pandas DataFrame where rows are sequence IDs and columns are positions with variation.
    """
    #create variation df
    variation_dict = {}

    #num_seqs = len(alignment)
    aln_length = alignment.get_alignment_length()
    
    #convert to numpy
    alignment_array = np.array([list(record.seq) for record in align])
    
    #get reference row for comparison
    ref_row = alignment_array[0]
    
    #find columns with variation
    is_variable = np.any(alignment_array != ref_row, axis=0)
    
    if nogaps:
        is_gap_free = ~np.any(alignment_array == '-', axis=0) #gap free columns 
        valid_columns = is_gap_free & is_variable
    else:
        valid_columns = is_variable
        
    #valid_positions = np.where(valid_columns)[0]
    valid_positions = np.flatnonzero(valid_columns) #faster?
    
    variation_df = pd.DataFrame(
        alignment_array[:, valid_columns],
        index=[record.id for record in alignment],
        columns=valid_positions
    )

    return variation_df