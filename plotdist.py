#!/usr/bin env python

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#FIXME:
from plot_functions.py import  plot_histogram()
############
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"


############
# Hist: All UP Proteomes
############
plot_histogram(data=df_up2,
                   x='protein_count',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='blue',
                   label='T', 
                   alpha=0.9,
                   ylog=True,
                   xlabel='Number of Proteins',
                   ylabel='Frequency',
                   extra_title_label='All UP Proteomes',
                   output_plot=basedir+"/Plots/HistAllUP_proteins.png")


plot_histogram(data=df_up2,
                   x='complete_combined_score',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='blue',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Busco',
                   ylabel='Frequency',
                   extra_title_label='All UP Proteomes',
                   output_plot=basedir+"/Plots/HistAllUP_BUSCO.png")

plot_histogram(data=df_up2,
                   x='checkm-completeness',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='blue',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='CheckM',
                   ylabel='Frequency',
                   extra_title_label='All UP Proteomes',
                   output_plot=basedir+"/Plots/HistAllUP_CheckM.png")


###############################################################################
############
# Hist: Redundant
############
redundant_df = df_up2[df_up2['is_redundant'].isin([-1, 1])]
#df_up2[df_up2['is_redundant'].isin([1, -1])]
redundant_df['is_redundant'].value_counts()

plot_histogram(data=redundant_df,
                   x='protein_count',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='grey',
                   label='T', 
                   alpha=0.9,
                   ylog=True,
                   xlabel='Number of Proteins',
                   ylabel='Frequency',
                   extra_title_label='Redundant UP Proteomes',
                   output_plot=basedir+"/Plots/HistRedundant_proteins.png")


plot_histogram(data=redundant_df,
                   x='complete_combined_score',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='grey',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Busco',
                   ylabel='Frequency',
                   extra_title_label='Redundant UP Proteomes',
                   output_plot=basedir+"/Plots/HistRedundant_BUSCO.png")

plot_histogram(data=redundant_df,
                   x='checkm-completeness',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='grey',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='CheckM',
                   ylabel='Frequency',
                   extra_title_label='Redundant UP Proteomes',
                   output_plot=basedir+"/Plots/HistRedundant_CheckM.png")


###############################################################################
############
# Hist: Excluded
############
excluded_df = df_up2[df_up2['is_excluded'] == 't']
excluded_df['is_excluded'].value_counts()

plot_histogram(data=excluded_df,
                   x='protein_count',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='red',
                   label='T', 
                   alpha=0.9,
                   ylog=True,
                   xlabel='Number of Proteins',
                   ylabel='Frequency',
                   extra_title_label='Excluded UP Proteomes',
                   output_plot=basedir+"/Plots/HistExcluded_proteins.png")


plot_histogram(data=excluded_df,
                   x='complete_combined_score',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='red',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Busco',
                   ylabel='Frequency',
                   extra_title_label='Excluded UP Proteomes',
                   output_plot=basedir+"/Plots/HistExcluded_BUSCO.png")

plot_histogram(data=excluded_df,
                   x='checkm-completeness',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='red',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='CheckM',
                   ylabel='Frequency',
                   extra_title_label='Excluded UP Proteomes',
                   output_plot=basedir+"/Plots/HistExcluded_CheckM.png")


###############################################################################
upids_cat= list(excluded_df['upid']) + list(redundant_df['upid'])
other_up = df_up2[~df_up2['upid'].isin(upids_cat)]
other_up.head()
other_up.shape

other_up2 = df_up2[~(df_up2['is_excluded'] == 't') & df_up2['is_redundant'].isin([-1, 1])]


plot_histogram(data=other_up,
                   x='protein_count',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='brown',
                   label='T', 
                   alpha=0.9,
                   ylog=True,
                   xlabel='Number of Proteins',
                   ylabel='Frequency',
                   extra_title_label='Other UP Proteomes',
                   output_plot=basedir+"/Plots/HistOtherlUP_proteins.png")


plot_histogram(data=other_up,
                   x='complete_combined_score',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='brown',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Busco',
                   ylabel='Frequency',
                   extra_title_label='Other UP Proteomes',
                   output_plot=basedir+"/Plots/HistOtherlUP_BUSCO.png")

plot_histogram(data=other_up,
                   x='checkm-completeness',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='brown',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='CheckM',
                   ylabel='Frequency',
                   extra_title_label='Other UP Proteomes',
                   output_plot=basedir+"/Plots/HistOtherlUP_CheckM.png")
