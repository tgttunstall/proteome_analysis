#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 14:21:24 2025

@author: tanu
"""
import os, sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#FIXME: the bins thing for the data
###############################################################################
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Input
# Data is coming from: plot_data_prep.py
# plot_data = atb_df_plot.copy()
#add_str="ATB"

plot_data = up_df_plot.copy()
add_str="UP"

print(len(plot_data.columns))
print(plot_data.shape)

outplot_checkm_scatter=basedir+"/Plots/Scatter_CM_CM2_" + add_str + ".png"
outplot_checkm_stats=basedir+"/Plots/Stats_CM_CM2_" + add_str + ".png"

outplot_CM_completeness_hist=basedir+"/Plots/Dist_completeness_CM_CM2_" + add_str + ".png"
outplot_CM_contamination_hist=basedir+"/Plots/Dist_contamination_CM_CM2_" + add_str + ".png"

###############################################################################
import pandas as pd
import scipy.stats as stats

# Assuming 'plot_data' is your DataFrame and 'Completeness_CM' is a column in your DataFrame
# Shapiro-Wilk Test
shapiro_test = stats.shapiro(plot_data['Completeness_CM'])
print(f"Shapiro-Wilk Test: Statistic={shapiro_test.statistic:.4f}, p-value={shapiro_test.pvalue:.4f}")

# D’Agostino’s K^2 Test
dagostino_test = stats.normaltest(plot_data['Completeness_CM'])
print(f"D’Agostino’s K^2 Test: Statistic={dagostino_test.statistic:.4f}, p-value={dagostino_test.pvalue:.4f}")


import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Set up the matplotlib figure
plt.figure(figsize=(18, 5))

# Histogram
plt.subplot(1, 3, 1)
sns.histplot(plot_data['Completeness_CM'], kde=True)
plt.title('Histogram of Completeness_CM')

# Q-Q plot
plt.subplot(1, 3, 2)
stats.probplot(plot_data['Completeness_CM'], dist="norm", plot=plt)
plt.title('Q-Q Plot of Completeness_CM')

# Box plot
plt.subplot(1, 3, 3)
sns.boxplot(x=plot_data['Completeness_CM'])
plt.title('Box Plot of Completeness_CM')

plt.tight_layout()
plt.show()
###############################################################################
#plot_data=plot_data[['Completeness_CM', 'Completeness_CM2', 'Contamination_CM', 'Contamination_CM2']]
#sns.pairplot(plot_data)
#sns.pairplot(plot_data, kind='kde')
#sns.pairplot(plot_data, corner=True)
###############################################################################

# colour
colour_cm = 'darkgoldenrod' #'blue' # cm
colour_cm2 = 'darkviolet' #darkorange #cm2

# font sizes
title_fontsize = 10
label_fontsize = 9
tick_fontsize = 9
legend_fontsize = 9

#============
#Hist plot OVERLAID: CM and CM2
#x=Contamination
#============
sns_plot_style='white' 

# alpha
cm_alpha = 0.8
cm2_alpha = 0.8

# data length
final_count = len(plot_data)

# hist bins
bin_method='rice'
nbins = np.histogram_bin_edges(plot_data['Contamination_CM'], bins=bin_method)
num_bins = len(nbins) - 1
bin_size = round(nbins[1] - nbins[0], 2) if num_bins > 1 else "Variable"

# Calculate mean and standard deviation for both datasets
mean_cm = plot_data['Contamination_CM'].mean()
std_cm = plot_data['Contamination_CM'].std()

mean_cm2 = plot_data['Contamination_CM2'].mean()
std_cm2 = plot_data['Contamination_CM2'].std()


sns.histplot(data=plot_data, 
             x='Contamination_CM', 
             bins=num_bins,
             #binwidth=50,
             color=colour_cm,
             label='CheckM', 
             alpha=cm_alpha, 
             #ax=axes[5]
             )
sns.histplot(data=plot_data,
             x='Contamination_CM2', 
             bins=num_bins,
             #binwidth=50,
             color=colour_cm2, 
             label='CheckM2', 
             alpha=cm2_alpha, 
             #ax=axes[5]
             )
plt.yscale('log')
plt.xlabel("Contamination (%)")
plt.ylabel("Frequency of Contamination")
plt.title(f'Distribution of Contamination in CM and CM2: {add_str} (n={n_data:,})\nBins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
plt.legend(fontsize=legend_fontsize)


# Add vertical lines for means
plt.axvline(mean_cm, color=colour_cm, linestyle='dashed', linewidth=2, label=f'CM Mean: {mean_cm:.2f}')
plt.axvline(mean_cm2, color=colour_cm2, linestyle='dashed', linewidth=2, label=f'CM2 Mean: {mean_cm2:.2f}')

# Add shading for ±1 SD
plt.fill_betweenx(y=[0, plt.gca().get_ylim()[1]], 
                  x1=mean_cm - std_cm, x2=mean_cm + std_cm, 
                  color=colour_cm, alpha=0.2, label=f'UP ±1 SD')

plt.fill_betweenx(y=[0, plt.gca().get_ylim()[1]], 
                  x1=mean_cm2 - std_cm2, x2=mean_cm2 + std_cm2, 
                  color=colour_cm2, alpha=0.2, label=f'ATB ±1 SD')

# additional summary stats:

# Existing code for calculation...
# Calculating median, interquartile range, and range for both datasets
median_cm = plot_data['Contamination_CM'].median()
iqr_cm = np.subtract(*np.percentile(plot_data['Contamination_CM'], [75, 25]))
range_cm = (plot_data['Contamination_CM'].max(), plot_data['Contamination_CM'].min())

median_cm2 = plot_data['Contamination_CM2'].median()
iqr_cm2 = np.subtract(*np.percentile(plot_data['Contamination_CM2'], [75, 25]))
range_cm2 = (plot_data['Contamination_CM2'].max(), plot_data['Contamination_CM2'].min())

final_count_cm = len(plot_data['Contamination_CM'])
final_count_cm2 = len(plot_data['Contamination_CM2'])
    
plt.tight_layout()

# Annotations for CM
plt.figtext(0.65, 0.65, 
            f"n = {final_count_cm:,}\nMean, SD = {mean_cm:.2f}, {std_cm:.2f}\nMedian, IQR = {median_cm:.2f}, {iqr_cm:.2f}\nRange = {range_cm[0]:.2f} - {range_cm[1]:.2f}", 
            fontsize=9, color=colour_cm, ha='left')


plt.figtext(0.65, 0.45, 
            f"n = {final_count_cm2:,}\nMean, SD = {mean_cm2:.2f}, {std_cm2:.2f}\nMedian, IQR = {median_cm2:.2f}, {iqr_cm2:.2f}\nRange = {range_cm2[0]:.2f} - {range_cm2[1]:.2f}", 
            fontsize=9, color=colour_cm2, ha='left')

# Save the figure
plt.savefig(outplot_CM_contamination_hist, format='png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()


plt.clf()
plt.close()
sns.set(style='white')  # Change the style to 'white'