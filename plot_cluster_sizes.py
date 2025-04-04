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
<<<<<<< HEAD
input_file  = basedir+"/Count_Labelled_Species_protein_cluster_atb.tsv"
=======
#input_file  = basedir+"/Count_Labelled_Species_protein_cluster.tsv"
>>>>>>> 18f27cb (added card1.py but need to move it)

# proteins
input_file  = basedir+"/Count_Proteins_Labelled_Species_protein_cluster_up.tsv"

##############################################################################

# Example usage:
<<<<<<< HEAD
input_df = pd.read_csv(input_file, sep='\t', names=['proteome_count'] )

# up+atb
#tot_proteomes = 146237
#tot_proteins = 307483396

# atb
tot_proteomes = 119701
tot_proteins = 252731905
=======
#input_df = pd.read_csv(input_file, sep='\t', names=['proteome_count'] )

input_df = pd.read_csv(input_file, sep='\t', names=['protein_count'] )


# UP + ATB
#tot_proteomes = 146237
#tot_proteins = 307483396

# UP
tot_proteomes = 26536
tot_proteins = 54751492
>>>>>>> 18f27cb (added card1.py but need to move it)


# X-axis: Number
plot_clustersize(df=input_df, 
#          plot_colname='proteome_count', 
          plot_colname='protein_count', 

          total_proteomes=tot_proteomes,
          total_proteins=tot_proteins,
          min_proteomes=None, 
          ylog=False, 
          x_proteomes='number',
<<<<<<< HEAD
          output_plot=basedir + "/Plots/Proteome_clustersizes_atb.png",
=======
#          output_plot=basedir + "/Plots/Proteome_clustersizes.png",
          output_plot=basedir + "/Plots/Protin_clustersizes_linear.png",
>>>>>>> 18f27cb (added card1.py but need to move it)
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
