#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch

#import random

homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
#up_file = basedir + "/up_pcounts.tsv"
up_file = basedir + "/spneumo_biosample_info_v1.out"
atb_file = basedir + "/atb_pcounts.tsv"

# Read the data
df_up_all = pd.read_csv(up_file, sep="\t")
#df_up_all = pd.read_csv(up_file, sep="\t", usecols = ['biosample', 'protein_count'])
df_atb_all = pd.read_csv(atb_file, sep="\t")

# sanity check: Drop duplicated, and data with protein count == 0
df_up = df_up_all[df_up_all['protein_count']> 0]
df_atb = df_atb_all[df_atb_all['protein_count'] > 0]

df_up2 = df_up.drop_duplicates(subset = 'biosample', keep = 'last')
df_up_dups = df_up[df_up.duplicated(subset='biosample')]

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
#categories = ["Only UP", "Both", "Only ATB"]
# counts = [len(only_up), len(both_sources), len(only_atb)]

# # Plot
# plt.figure(figsize=(8, 7))
# plt.bar(categories, counts, color=["blue", "green", "red"], alpha=0.7)

# # Labels & Title
# plt.xlabel("Data Sources")
# plt.ylabel("Number of Proteomes")
# plt.title("Comparison of Proteome counts by Biosample in UP and ATB")
# plt.grid(axis="y", linestyle="--", alpha=0.5)

# # Annotate counts
# for i, v in enumerate(counts):
#     plt.text(i, v + max(counts) * 0.02, str(v), ha="center", fontsize=14)

# plt.show()
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
venn.get_patch_by_id("01").set_color("red")
venn.get_patch_by_id("11").set_color("green")

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
            color=["blue", "green", "red"],
            alpha=0.2, 
            #color=["lightgrey", "lightgrey", "lightgrey"],
            width=bar_width)
# Overlay bars (colored with denser hatching for exclusive)
axes[1].bar(categories, 
            subset_counts, 
            color=["blue", "green", "red"], 
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

###############################################################################
#=========================
# Histogram: Protein count distribution
#=========================
# Extract counts
up_proteins = df_up2["protein_count"]
atb_proteins = df_atb2["protein_count"]

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

# --- PLOT: histogram with bin width specified  ---

# Define bins (e.g., group by protein count ranges)
bins = np.arange(1, max(up_proteins.max(), atb_proteins.max()) + 100, 100)  # Adjust bin size as needed

# Count occurrences per bin
up_counts, _ = np.histogram(up_proteins, bins)
atb_counts, _ = np.histogram(atb_proteins, bins)

# Bar plot
plt.figure(figsize=(8, 6))
plt.bar(bins[:-1], 
        up_counts, 
        width=100, 
        align="edge", 
        alpha=0.6, 
        color="blue", 
        label="UP", 
        edgecolor="black")

plt.bar(bins[:-1], 
        atb_counts, 
        width=100, 
        align="edge", 
        alpha=0.6, 
        color="red", 
        label="ATB", 
        edgecolor="black")

# Apply log scale to the y-axis
plt.yscale("log")
#plt.xscale("log")

# Labels & Title
plt.xlabel("Number of Proteins")
plt.ylabel("Number of Proteomes")
plt.title("Protein Count Distribution Across Proteomes")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.show()
