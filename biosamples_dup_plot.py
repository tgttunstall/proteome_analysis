#!/bin/env python

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

###############################################################################
#FIXME:
from plot_functions.py import  plot_protein_busco()
############
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"
outdir_biosample_plots = basedir+"/Plots/biosamples"
#########
#Data for plotting
#df_pivot
def chunk_dataframe(df, chunk_size):
    #chunk_size = len(df) // chunks
    #return [df[i:i+chunk_size] for i in range(0, len(df), chunk_size)]
    return [df[i:i+chunk_size].reset_index(drop=True) for i in range(0, len(df), chunk_size)]


chunked_dfs = chunk_dataframe(df=df_pivot, chunk_size=10,reset_index=True)

for i, chunk in enumerate(chunked_dfs):
    print(f"i: {i}, chunk:\n {chunk}")
    # Plot or process each chunk
    f, a, a2 = plot_protein_busco(chunk, output_plot=outdir_biosample_plots +"/" +f"Chunk {i+1}.png")



#########
fig, ax1, ax2 = plot_protein_busco(
    #df=df_pivot.iloc[0:10,],
    df=c2,
    #df=chunk,
    spacing_factor=0.85,
    bar_width= 0.15,
    gap=0.18,
    title_fontsize=14, 
    label_fontsize=14, 
    tick_fontsize=12, 
    legend_fontsize=14, 
    fig_height = 11,
    plot_title="Duplicate Biosamples: Protein count and BUSCO",
    #output_plot=outdir_biosample_plots + "/T1.png"
    )
plt.show()
