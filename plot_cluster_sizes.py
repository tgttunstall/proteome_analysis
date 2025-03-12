#!/bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

#FIXME:
#from plot_functions.py import  plot_clustersize()
###############################################################################
# Read data
homedir = os.path.expanduser("~")
basedir =  homedir + "/Documents/arise/spneumo_dataset"
#basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

#input_file  = basedir+"/DEL/clusterizes_proteomes2"
input_file  = basedir+"/Count_Labelled_Species_protein_cluster_atb.tsv"


##############################################################################

# Example usage:
input_df = pd.read_csv(input_file, sep='\t', names=['proteome_count'] )

# up+atb
#tot_proteomes = 146237
#tot_proteins = 307483396

# atb
tot_proteomes = 119701
tot_proteins = 252731905


# X-axis: Number
plot_clustersize(df=input_df, 
          plot_colname='proteome_count', 
          total_proteomes=tot_proteomes,
          total_proteins=tot_proteins,
          min_proteomes=None, 
          ylog=True, 
          x_proteomes='number',
          output_plot=basedir + "/Plots/Proteome_clustersizes_atb.png",
          show_stats=True)

plot_clustersize(df=input_df, 
          plot_colname='proteome_count', 
          total_proteomes=tot_proteomes,
          total_proteins=tot_proteins,
          min_proteomes=None, 
          ylog=True, 
          x_proteomes='number',
          output_plot=basedir + "/Plots/Proteome_clustersizes_filtered_atb.png",
          show_stats=True)


# X-axis: percentage
plot_clustersize(df = input_df, 
          plot_colname = 'proteome_count', 
          total_proteomes = tot_proteomes,
          total_proteins = tot_proteins,
          min_proteomes = None, 
          ylog = True, 
          x_proteomes = 'percent',
          output_plot = basedir + "/Plots/Proteome_clustersizes_percent_atb.png",
          show_stats = True)

plot_clustersize(df = input_df, 
          plot_colname = 'proteome_count', 
          total_proteomes = tot_proteomes,
          total_proteins = tot_proteins,
          min_proteomes = 1, 
          ylog = True, 
          x_proteomes = 'percent', 
          output_plot = basedir + "/Plots/Proteome_clustersizes_percent_filtered_atb.png",
          show_stats = True)
