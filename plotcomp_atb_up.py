#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import seaborn as sns


#import random
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
#up_file = basedir + "/up_pcounts.tsv"
#up_file = basedir + "/spneumo_biosample_info_v1.out"
up_file = basedir + "/spneumo_biosample_info_processed.out"
atb_file = basedir + "/atb_pcounts.tsv"

# Read the data
df_up_all = pd.read_csv(up_file, sep="\t")
#df_up_all = pd.read_csv(up_file, sep="\t", usecols = ['biosample', 'protein_count'])
df_atb_all = pd.read_csv(atb_file, sep="\t")

# sanity check: Drop duplicated, and data with protein count == 0
df_up = df_up_all[df_up_all['protein_count']> 0]
df_atb = df_atb_all[df_atb_all['protein_count'] > 0]

df_up2 = df_up.drop_duplicates(subset = ['biosample'], keep = 'last')
df_up_dups = df_up[df_up.duplicated(subset= ['biosample', 'upid'])]

df_up_dups.to_csv(basedir + "/up_dups.tsv", sep = "\t")

df_atb_dups = df_atb[df_atb.duplicated(subset = 'biosample')]
df_atb2 = df_atb.drop_duplicates(subset=['biosample'], keep = 'last')

# Merge data on 'biosample' to see coverage
#df_merged = pd.merge(df_up, df_atb, on="biosample", how="outer", suffixes=("_up", "_atb"))
df_merged = pd.merge(df_up2, df_atb2, on="biosample", how="outer", suffixes=("_up", "_atb"))

# Replace NaN with 0 (if a biosample is missing in one dataset)
df_merged.fillna(0, inplace=True)

# Count different categories
only_up = df_merged[(df_merged["protein_count_up"] > 0) & (df_merged["protein_count_atb"] == 0)]
only_atb = df_merged[(df_merged["protein_count_atb"] > 0) & (df_merged["protein_count_up"] == 0)]
both_sources = df_merged[(df_merged["protein_count_up"] > 0) & (df_merged["protein_count_atb"] > 0)]

###############################################################################
#===========
# Bar plot
#===========
# --- PLOT: matplotlib ---
categories = ["Only UP", "Both", "Only ATB"]
counts = [len(only_up), len(both_sources), len(only_atb)]

plt.figure(figsize=(8, 7))
plt.bar(categories, counts, color=["blue", "peru", "green"], alpha=0.7)

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Proteomes")
plt.title("Comparison of Proteome counts by Biosample in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts
for i, v in enumerate(counts):
    plt.text(i, v + max(counts) * 0.02, str(v), ha="center", fontsize=14)

plt.show()

# --- PLOT: Seaborn ---
categories = ["Only UP", "Both", "Only ATB"]
counts = [len(only_up), len(both_sources), len(only_atb)]

data = pd.DataFrame({
    'Data Sources': categories,
    'Number of Proteomes': counts
    })

# Plotting using Seaborn
plt.figure(figsize=(8, 7))
bar_plot = sns.barplot(x='Data Sources', y='Number of Proteomes', data=data, alpha=0.7, palette=["blue", "peru", "green"])

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Proteomes")
plt.title("Comparison of Proteome Counts by Biosample in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts
for index, row in data.iterrows():
    bar_plot.text(index, row['Number of Proteomes'] + max(counts) * 0.01, str(row['Number of Proteomes']), color='black', ha="center", fontsize=14)

plt.show()

###############################################################################
#===============
# Venn and Stacked bar plot: side by side
#===============
# Define set sizes
only_up_count = len(only_up)  # Unique to UP
only_atb_count = len(only_atb)  # Unique to ATB
both_count = len(both_sources)  # In both UP and ATB

# Define bar plot data
categories = ["UP", "UP+ATB", "ATB"]
total_counts = [only_up_count + both_count, 
                only_up_count + only_atb_count + both_count,
                only_atb_count + both_count]
subset_counts = [only_up_count, both_count, only_atb_count]

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# --- Venn Diagram ---
venn = venn2(subsets=(only_up_count, only_atb_count, both_count), set_labels=("UP", "ATB"), ax=axes[0])

# Customize Venn colors
venn.get_patch_by_id("10").set_color("blue")
venn.get_patch_by_id("01").set_color("green")
venn.get_patch_by_id("11").set_color("peru")

# Annotate numbers
venn.get_label_by_id("10").set_text(f"{only_up_count}")
venn.get_label_by_id("01").set_text(f"{only_atb_count}")
venn.get_label_by_id("11").set_text(f"{both_count}")

axes[0].set_title("Comparing Proteomes/Biosamples in UP and ATB")

# --- Stacked Bar Plot ---
bar_width = 0.6

# Base bars (colored)
axes[1].bar(categories, 
            total_counts, 
            color=["blue", "peru", "green"],
            alpha=0.2, 
            #color=["lightgrey", "lightgrey", "lightgrey"],
            width=bar_width)
# Overlay bars (colored with denser hatching for exclusive)
axes[1].bar(categories, 
            subset_counts, 
            color=["blue", "peru", "green"], 
            alpha=0.8, 
            width=bar_width, 
            hatch='////')  # Denser hatching

# Labels & Title
axes[1].set_xlabel("Data Sources")
axes[1].set_ylabel("Number of Proteomes")
#axes[1].set_title("Stacked Bar Plot of Proteome Counts")
axes[1].grid(axis="y", 
             linestyle="--", 
             alpha=0.5)

# Annotate counts on bars
for i, (total, subset) in enumerate(zip(total_counts, subset_counts)):
    axes[1].text(i, 
                 total + max(total_counts) * 0.01, 
                 str(total), 
                 ha="center", 
                 fontsize=11, 
                 color="black")
    axes[1].text(i, subset / 2, str(subset), 
                 ha="center", 
                 fontsize=11, 
                 color="white", 
                 fontweight="bold")  # Centered inside subset bar

#axes[1].legend(["Total", "Exclusive"])
# Custom legend
legend_elements = [
    Patch(facecolor='lightgrey', label='Total'),
    Patch(facecolor='lightgrey', hatch='////', label='Exclusive')  # Denser hatch
]
axes[1].legend(handles=legend_elements)

# Show both plots
plt.tight_layout()
plt.show()

# 25961 up

###############################################################################
#=========================
# Histogram: Protein count distribution
#=========================
# Extract counts
# up_proteins = df_up2["protein_count"] # dups removed based on biosample for accurate atb count overlap
up_proteins = df_up["protein_count"]
atb_proteins = df_atb2["protein_count"]

df_up["protein_count"].describe()
df_atb2["protein_count"].describe()

# --- PLOT: histogram with number of bins specified --- #not always accurate!
# plt.figure(figsize=(8, 6))
# plt.hist(up_proteins, bins=50, alpha=0.6, color="blue", label="UP", edgecolor="black")
# plt.hist(atb_proteins, bins=50, alpha=0.6, color="red", label="ATB", edgecolor="black")

# # Apply log scale to the y-axis
# plt.yscale("log")
# #plt.xscale("log")

# # Labels & Title
# plt.xlabel("Number of Proteins")
# plt.ylabel("Number of Proteomes")
# plt.title("Protein Count Distribution Across Proteomes")
# plt.legend()
# plt.grid(axis="y", linestyle="--", alpha=0.5)
# plt.show()

# --- PLOT: histogram with bin width specified  (Matplotlib) ---
# Define bins (e.g., group by protein count ranges)
bin_width = 100 
bins = np.arange(1, max(up_proteins.max(), atb_proteins.max()) + bin_width, bin_width)  # Adjust bin size as needed

# Count occurrences per bin
up_counts, _ = np.histogram(up_proteins, bins)
atb_counts, _ = np.histogram(atb_proteins, bins)

num_bins_up = (up_counts>0).sum()
num_bins_atb = (atb_counts>0).sum()


###############################################################################
#=========================
# Bar plot (Matplotlib)
#=========================
plt.figure(figsize=(8, 6))
plt.bar(bins[:-1], 
        up_counts, 
        width=100, 
        align="edge", 
        alpha=0.6, 
        color="blue", 
        #label="UP",
        label=f'UP (Bin size: {bin_width}, Bin count: {num_bins_atb})',
        edgecolor="black")

plt.bar(bins[:-1], 
        atb_counts, 
        width=100, 
        align="edge", 
        alpha=0.6, 
        color="green", 
        #label="ATB", 
        label=f'ATB (Bin size: {bin_width}, Bin count: {num_bins_atb})',
        edgecolor="black")

# Apply log scale to the y-axis
plt.yscale("log")
#plt.xscale("log")

# Labels & Title
plt.xlabel("Number of Proteins")
plt.ylabel("Number of Proteomes")
plt.title("Protein Count Distribution: UP vs ATB Proteomes")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.show()
###############################################################################
# --- PLOT: histogram with bin width specified  (Seaborn) ---
# Define larger font sizes for better readability
title_fontsize = 15
label_fontsize = 13
tick_fontsize = 12
legend_fontsize = 11

#num_bins_up = (up_counts>0).sum()
#num_bins_atb = (atb_counts>0).sum()


df_up["protein_count"].describe()
df_atb2["protein_count"].describe()



# Setup the plotting
plt.figure(figsize=(8, 6))

# Histogram for All UP data
sns.histplot(#data=df_up2, 
             data=df_up,
             x='protein_count', 
             bins=bins, 
             color='blue',
             #label=f'UP (n={len(df_up2)}, Bin size: {bin_width}, Bin count: {num_bins_up})',
             label=f'UP (n={len(df_up)}, Bin count: {num_bins_up})',
             alpha=0.5, 
             #element="step", 
             edgecolor='black', 
             linewidth=1)

# Histogram for ATB data
sns.histplot(data=df_atb2, 
             x='protein_count', 
             bins=bins, 
             color='green',
             #label=f'ATB (n={len(df_atb2)}, Bin size: {bin_width}, Bin count: {num_bins_atb})',
             label=f'ATB (n={len(df_atb2)}, Bin count: {num_bins_atb})',

             alpha=0.5, 
             #element="step",
             edgecolor = 'k',
             linewidth=1)

# Adjusting plot aesthetics
plt.yscale('log')
plt.xlabel(f'Number of Proteins (Bin size: {bin_width})', fontsize=label_fontsize)
plt.ylabel('Number of Proteomes', fontsize=label_fontsize)
plt.title('Protein Count Distribution: UP vs ATB Proteomes', fontsize=title_fontsize)
plt.tick_params(axis='both', labelsize=tick_fontsize)
plt.legend(title='', fontsize=legend_fontsize)
plt.show()

################################################################################
# --- PLOT: histogram with bin width specified  (All UP data with categories) (Seaborn)
# Plotting setup
title_fontsize = 16
label_fontsize = 14
tick_fontsize = 12
legend_fontsize = 12

general_alpha=0.2
redundant_alpha=0.3
excluded_alpha=0.2
reference_alpha=0.8
rep_alpha=0.8

fig, axes = plt.subplots(2, 1, figsize=(12, 16), sharex=True)  # 2 rows, 1 column, sharing x-axis

# Top plot for df_up2 with different categories
sns.histplot(data=df_up2, 
             x='protein_count', 
             bins=bins, 
             color='blue', 
             label='General UP Proteomes', 
             alpha=general_alpha, 
             ax=axes[0])

sns.histplot(data=df_up2[df_up2['is_redundant'].isin([1, -1])], 
             x='protein_count', 
             bins=bins, 
             color='orange', 
             label='Redundant', 
             alpha=redundant_alpha, 
             ax=axes[0])

sns.histplot(data=df_up2[df_up2['is_excluded'] == 't'], 
             x='protein_count', 
             bins=bins, 
             color='red', 
             label='Excluded', 
             alpha=excluded_alpha, 
             ax=axes[0])

sns.histplot(data=df_up2[df_up2['is_reference'] == 't'], 
             x='protein_count', 
             bins=bins, 
             color='magenta', 
             label='Reference', 
             alpha=reference_alpha, 
             hatch="--", 
             ax=axes[0])

sns.histplot(data=df_up2[df_up2['is_representative'] == 't'], 
             x='protein_count', 
             bins=bins, 
             color='cyan', 
             label='Representative', 
             alpha=rep_alpha, 
             hatch="////", 
             ax=axes[0])

axes[0].set_yscale('log')
axes[0].set_xlabel('Number of Proteins', fontsize=label_fontsize)
axes[0].set_ylabel('Number of Proteomes', fontsize=label_fontsize)
axes[0].set_title('UP Proteome Distribution by Category', fontsize=title_fontsize)
axes[0].tick_params(axis='both', labelsize=tick_fontsize)
axes[0].legend(fontsize=legend_fontsize)

# Bottom plot for df_atb2 with general counts
sns.histplot(data=df_atb2, x='protein_count', bins=bins, color='green', label='ATB Proteomes', alpha=0.5, ax=axes[1])

axes[1].set_yscale('log')
axes[1].set_xlabel('Number of Proteins', fontsize=label_fontsize)
axes[1].set_ylabel('Number of Proteomes', fontsize=label_fontsize)
axes[1].set_title('ATB Proteome Count', fontsize=title_fontsize)
axes[1].tick_params(axis='both', labelsize=tick_fontsize)
axes[1].legend(fontsize=legend_fontsize)

plt.tight_layout()
plt.show()
##############################################################################

# --- PLOT: histogram with bin width specified  (Filtered UP data) (Seaborn)
# Filter out 'is_excluded' marked as 't' from df_up2
filtered_df_up2 = df_up2[df_up2['is_excluded'] != 't']

# Define bins for the histograms, considering both datasets
max_count = max(filtered_df_up2['protein_count'].max(), df_atb2['protein_count'].max())
bins = np.arange(1, max_count + bin_width, bin_width)
#num_bins_up = len(bins) - 1
#num_bins_atb = len(bins) - 1

up_counts_filtered, _ = np.histogram(filtered_df_up2['protein_count'], bins)
atb_counts, _ = np.histogram(df_atb2['protein_count'], bins)
#len(bins) - (up_counts_filtered==0).sum() -1 
#len(bins) - (atb_counts==0).sum() -1 

num_bins_up = (up_counts_filtered>0).sum()
num_bins_atb = (atb_counts>0).sum()

# Setup the plotting
plt.figure(figsize=(12, 8))

# Histogram for filtered UP data
sns.histplot(data=filtered_df_up2, 
             x='protein_count', 
             bins=bins, 
             color='blue',
             label=f'Filtered UP (Bin size: {bin_width}, Bin count: {num_bins_up})',
             alpha=0.5, 
             #element="step", 
             edgecolor='black', 
             linewidth=1)

# Histogram for ATB data
sns.histplot(data=df_atb2, 
             x='protein_count', 
             bins=bins, 
             color='green',
             label=f'ATB (Bin size: {bin_width}, Bin count: {num_bins_atb})',
             alpha=0.5, 
             #element="step",
             edgecolor = 'k',
             linewidth=1)

# Adjusting plot aesthetics
plt.yscale('log')
plt.xlabel('Number of Proteins', fontsize=14)
plt.ylabel('Number of Proteomes', fontsize=14)
plt.title('Protein Count Distribution: Filtered UP vs ATB', fontsize=title_fontsize)
plt.tick_params(axis='both', labelsize=tick_fontsize)

plt.legend(title='', fontsize = legend_fontsize)
plt.show()
