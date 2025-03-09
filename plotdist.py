#!/usr/bin env python

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#FIXME:
from plot_functions.py import  plot_histogram()

############

plot_histogram(data=df_up2,
                   x='protein_count',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='blue',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Number of Proteins',
                   ylabel='Frequency',
                   output_plot=None)

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
                   output_plot=None)

plot_histogram(data=df_up2,
                   x='checkm-completeness',
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='blue',
                   label='T', 
                   alpha=0.5,
                   ylog=True,
                   xlabel='Busco',
                   ylabel='Frequency',
                   output_plot=None)
