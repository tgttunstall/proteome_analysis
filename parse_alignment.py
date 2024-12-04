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

ALIGN_FORMAT = "fasta"

homedir = os.path.expanduser("~")

# time this with diff formats: clustal, fasta, etc
#aln_fname = "/home/tunstall/git/arise_proteome/11576.aln"
#aln_fname = "/home/pub/Work/data_arise_proteome/testc2/.aln"
#aln_fname = "/home/tunstall/git/arise_proteome/829.aln"

aln_fname = "/home/pub/Work/data_arise_proteome/testC/results_testC_v3/20613.aln"

align = AlignIO.read(aln_fname, ALIGN_FORMAT)
align[:,0]


variation_dict = {}

aln_length = len(align[0])
for pos in range(0, aln_length):
    print(f"position:{pos}")
    column = align[:, pos]
    print(f"column:{column}")

    prev_res = column[0]
    if prev_res == "-":
        continue
                        
    for res in column[1:]:
        if res == "-":
            break
        if res != prev_res:
            print(f"{pos} {column}")
            variation_dict[pos] = list(column)
            break

#dict to df
#add column with all the ids
variation_df = pd.DataFrame.from_dict(variation_dict, orient = "columns")
variation_df['id'] = [align[i].id for i in range(len(variation_df))]

# make id the first column:Bho
variation_df = variation_df[['id'] + [col for col in variation_df.columns if col != 'id']]
#ids = [align[i].id for i in range(len(variation_df))]

#or if you want that as index
#variation_df = pd.DataFrame.from_dict(variation_dict, orient = "columns")
#variation_df.index = [align[i].id for i in range(len(align))]

variation_df
#
#variation_df.to_csv("varcsv, ")

#######################
#that prints the columns for variations
#we need to transpose it and then add the headers for the columns captured (which are actially the rows)


#
#id               221 223
#35497:ENSCP....  K M
#123:ENSCP...     K N
#2323:ENSCP..     L P