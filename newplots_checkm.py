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

###############################################################################
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Input
# Data is coming from: plot_data_prep.py
#plot_data = atb_df_plot.copy()
#add_str="ATB"

plot_data = up_df_plot.copy()
add_str="UP"

print(len(plot_data.columns))
print(plot_data.shape)

outplot_checkm_scatter=basedir+"/Plots/Scatter_CM_CM2_" + add_str + ".png"
outplot_checkm_stats=basedir+"/Plots/Stats_CM_CM2_" + add_str + ".png"

outplot_CM_completeness_hist=basedir+"/Plots/Dist_completeness_CM_CM2_" + add_str + ".png"
outplot_CM_contamination_hist=basedir+"/Plots/Dist_contamination_CM_CM2_" + add_str + ".png"

outplot_checkm_busco=basedir+"/Plots/Scatter_BUSCO_CM_CM2_" + add_str + ".png"

###############################################################################
#plot_data=plot_data[['Completeness_CM', 'Completeness_CM2', 'Contamination_CM', 'Contamination_CM2']]
#sns.pairplot(plot_data)
#sns.pairplot(plot_data, kind='kde')
#sns.pairplot(plot_data, corner=True)
###############################################################################

# colour
colour_cm = 'darkgoldenrod' # cm
colour_cm2 = 'darkviolet' #cm2

up_col = 'blue'
atb_col = 'darkorange'

# font sizes
title_fontsize = 10
label_fontsize = 9
tick_fontsize = 9
legend_fontsize = 9

#===============
#Scatter plot : checkM
#x=Completeness
#y=Contamination
#===============
#------------
# matplotlib
#------------
#ax1 = plot_data.plot(kind='scatter', y='Completeness_CM', x='Contamination_CM', color=colour_cm, label='CM')
#ax2 = plot_data.plot(kind='scatter', y='Completeness_CM2', x='Contamination_CM2', color=colour_cm2', label='CM2', ax=ax1)

# Specify x-axis and y-axis labels
#ax1.set_ylabel('Completeness')
#ax1.set_xlabel('Contamination')

#------------
# Seaborn
#------------
sns_plot_style='white'  # whitegrid, dark, darkgrid, ticks

legend_pos='upper left'
n_data=len(plot_data)

plt.figure(figsize=(10, 6))
sns.scatterplot(data=plot_data, x='Completeness_CM', y='Contamination_CM', color=colour_cm, label='CheckM')
sns.scatterplot(data=plot_data, x='Completeness_CM2',y='Contamination_CM2', color=colour_cm2, label='CheckM2')

# Adding labels and title
plt.xlabel('Completeness (%)')
plt.ylabel('Contamination (%)')
plt.title(f'Comparison of CheckM and CheckM2: {add_str} (n={n_data:,})')

# Display the legend
plt.legend(loc=legend_pos)

# Add stats
#plt.figtext(0.15, 0.55, f"n = {final_count_cm}\nMean, SD = {mean_cm:.1f}, {std_cm:.1f}\nMedian, IQR = {median_cm:.1f}, {iqr_cm:.1f}\nRange = {range_cm[0]:.1f} - {range_cm[1]:.1f}", 
#            fontsize=9, color=colour_cm, ha='left')

#plt.figtext(0.15, 0.40, f"n = {final_count_cm2}\nMean, SD = {mean_cm2:.1f}, {std_cm2:.1f}\nMedian, IQR = {median_cm2:.1f}, {iqr_cm2:.1f}\nRange = {range_cm2[0]:.1f} - {range_cm2[1]:.1f}", 
#            fontsize=9, color=colour_cm2, ha='left')

# Save the figure
plt.savefig(outplot_checkm_scatter, format='png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
###############################################################################
#===============
#Stats plot: To go with scatter plot
#Completeness
#Contamination
#===============
# Calculate statistics for both datasets
comp_stats_cm = {
    'mean': plot_data['Completeness_CM'].mean(),
    'std': plot_data['Completeness_CM'].std(),
    'median': plot_data['Completeness_CM'].median(),
    'iqr': plot_data['Completeness_CM'].quantile(0.75) - plot_data['Completeness_CM'].quantile(0.25),
    'range': (plot_data['Completeness_CM'].max(), plot_data['Completeness_CM'].min())
}

comp_stats_cm2 = {
    'mean': plot_data['Completeness_CM2'].mean(),
    'std': plot_data['Completeness_CM2'].std(),
    'median': plot_data['Completeness_CM2'].median(),
    'iqr': plot_data['Completeness_CM2'].quantile(0.75) - plot_data['Completeness_CM2'].quantile(0.25),
    'range': (plot_data['Completeness_CM2'].max(), plot_data['Completeness_CM2'].min())
}

# Calculate statistics for Contamination_CM
contam_stats_cm = {
    'mean': plot_data['Contamination_CM'].mean(),
    'std': plot_data['Contamination_CM'].std(),
    'median': plot_data['Contamination_CM'].median(),
    'iqr': plot_data['Contamination_CM'].quantile(0.75) - plot_data['Contamination_CM'].quantile(0.25),
    'range': (plot_data['Contamination_CM'].max(), plot_data['Contamination_CM'].min())
}

# Calculate statistics for Contamination_CM2
contam_stats_cm2 = {
    'mean': plot_data['Contamination_CM2'].mean(),
    'std': plot_data['Contamination_CM2'].std(),
    'median': plot_data['Contamination_CM2'].median(),
    'iqr': plot_data['Contamination_CM2'].quantile(0.75) - plot_data['Contamination_CM2'].quantile(0.25),
    'range': (plot_data['Contamination_CM2'].max(), plot_data['Contamination_CM2'].min())
}

# Data for tables
data_completeness = [
    ["Mean, SD", f"{comp_stats_cm['mean']:.2f}, {comp_stats_cm['std']:.2f}", f"{comp_stats_cm2['mean']:.2f}, {comp_stats_cm2['std']:.2f}"],
    ["Median, IQR", f"{comp_stats_cm['median']:.2f}, {comp_stats_cm['iqr']:.2f}", f"{comp_stats_cm2['median']:.2f}, {comp_stats_cm2['iqr']:.2f}"],
    ["Range", f"{comp_stats_cm['range'][0]:.2f} - {comp_stats_cm['range'][1]:.2f}", f"{comp_stats_cm2['range'][0]:.2f} - {comp_stats_cm2['range'][1]:.2f}"]
]

data_contamination = [
    ["Mean, SD", f"{contam_stats_cm['mean']:.2f}, {contam_stats_cm['std']:.2f}", f"{contam_stats_cm2['mean']:.2f}, {contam_stats_cm2['std']:.2f}"],
    ["Median, IQR", f"{contam_stats_cm['median']:.2f}, {contam_stats_cm['iqr']:.2f}", f"{contam_stats_cm2['median']:.2f}, {contam_stats_cm2['iqr']:.2f}"],
    ["Range", f"{contam_stats_cm['range'][0]:.2f} - {contam_stats_cm['range'][1]:.2f}", f"{contam_stats_cm2['range'][0]:.2f} - {contam_stats_cm2['range'][1]:.2f}"]
]

# Create figure
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4,3))
ax1.axis('off')
ax2.axis('off')

# Adjust subplot spacing
fig.subplots_adjust(hspace=0.00001)  # Adjust the horizontal space between the tables

# Completeness Table
table1 = ax1.table(cellText=data_completeness,
                  colLabels=['Completeness', 'CheckM', 'CheckM2'],
                  loc='center',
                  cellLoc='center',
                  colColours=["dimgray", colour_cm, colour_cm2])  # CM and CM2 color scheme
table1.auto_set_font_size(False)
table1.set_fontsize(10)
table1.scale(1, 1.2)  # Adjust table scaling to fit

# Contamination Table
table2 = ax2.table(cellText=data_contamination,
                  colLabels=['Contamination', 'CheckM', 'CheckM2'],
                  loc='center',
                  cellLoc='center',
                  colColours=["lightgray", colour_cm, colour_cm2])  # Consistent with CM and CM2 colors
table2.auto_set_font_size(False)
table2.set_fontsize(10)
table2.scale(1, 1.2)  # Adjust table scaling to fit

plt.suptitle(f'Summary Stats for CM and CM2: {add_str} (n={n_data:,})', fontsize=title_fontsize)

# Save the figure
plt.savefig(outplot_checkm_stats, format='png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
###############################################################################
#===============
#Scatter plot : BUSCO vs CheckM
#x=BUSCO
#y=CheckM completeness and CheckM2 completeness
#===============
# Correlation
pearson_corr_cm = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM'], method='pearson')
spearman_corr_cm  = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM'], method='spearman')
kendall_corr_cm = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM'], method='kendall')

pearson_corr_cm2 = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM2'], method='pearson')
spearman_corr_cm2 = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM2'], method='spearman')
kendall_corr_cm2 = plot_data['complete_combined_score'].corr(plot_data['Completeness_CM2'], method='kendall')
#------------
# Seaborn
#------------
sns_plot_style='white' 

legend_pos='upper left'
n_data=len(plot_data)

plt.figure(figsize=(6, 6))
sns.scatterplot(data=plot_data, x='complete_combined_score', y='Completeness_CM', color=colour_cm, label='CM') #label=add_str+'_CM'
sns.scatterplot(data=plot_data, x='complete_combined_score',y='Completeness_CM2', color=colour_cm2, label='CM2')

# Adding labels and title
plt.xlabel('BUSCO (%)')
plt.ylabel('CheckM_completeness (%)')
plt.title(f'Comparison of BUSCO and CheckM: {add_str} (n={n_data:,})')

# Display the legend
plt.legend(loc=legend_pos)

# Annotate the plot with correlation coefficients
correlation_text_cm = f"Pearson: {pearson_corr_cm:.2f}\nSpearman: {spearman_corr_cm:.2f}\nKendall: {kendall_corr_cm:.2f}"
correlation_text_cm2 = f"Pearson: {pearson_corr_cm2:.2f}\nSpearman: {spearman_corr_cm2:.2f}\nKendall: {kendall_corr_cm2:.2f}"

plt.figtext(0.15, 0.70, correlation_text_cm, fontsize=10, color=colour_cm, ha='left')
plt.figtext(0.15, 0.59, correlation_text_cm2, fontsize=10, color=colour_cm2, ha='left')

# Save the figure
plt.savefig(outplot_checkm_busco, format='png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
###############################################################################
#RESET if needed something else
plt.clf()
plt.close()
sns.set(style='white')  # Change the style to 'white'
###############################################################################