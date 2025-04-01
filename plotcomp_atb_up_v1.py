#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import seaborn as sns
#import itertools


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

# ATB
df_atb_dups = df_atb[df_atb.duplicated(subset = 'biosample')]
print(list(df_atb_dups['biosample']))
#c = ['SAMN11121455','SAMEA115795954', 'SAMEA4024873D','SAMN03765066'] #sanity check
df_atb_dups= df_atb[df_atb['biosample'].isin(list(df_atb_dups['biosample']))]
print(df_atb_dups)
df_atb2 = df_atb.drop_duplicates(subset=['biosample'], keep = 'first')

# UP
df_up_biosample_dups = df_up[df_up.duplicated(subset = ['biosample'], keep = False)]
n_up_biosample_dups = df_up_biosample_dups['biosample'].nunique()
print(f"\nLength of up df with duplicate biosamples: {len(df_up_biosample_dups)}")
print(f"Length of unique UP duplicated biosamples: {n_up_biosample_dups}")
#df_up_biosample_dups.to_csv(basedir + "/up_biosample_dups.tsv", sep = "\t")

#------- CHECK
df_up2 = df_up.drop_duplicates(subset = ['biosample'], keep = 'last')
#df_up2 = df_up.copy()

# Merge data on 'biosample' to see coverage
#df_merged = pd.merge(df_up, df_atb, on="biosample", how="outer", suffixes=("_up", "_atb"))
#df_merged = pd.merge(df_up, df_atb2, on="biosample", how="outer", suffixes=("_up", "_atb"))
df_merged = pd.merge(df_up2, df_atb2, on="biosample", how="outer", suffixes=("_up", "_atb"))
# Replace NaN with 0 (if a biosample is missing in one dataset)
df_merged.fillna(0, inplace=True)

#------- CHECK
#df_merged1 = pd.merge(df_up, df_atb2, on="biosample", how="outer", suffixes=("_up", "_atb"))
#df_merged1.fillna(0, inplace=True)


# Count different categories
only_up = df_merged[(df_merged["protein_count_up"] > 0) & (df_merged["protein_count_atb"] == 0)]
only_atb = df_merged[(df_merged["protein_count_atb"] > 0) & (df_merged["protein_count_up"] == 0)]
both_sources = df_merged[(df_merged["protein_count_up"] > 0) & (df_merged["protein_count_atb"] > 0)]
both_sources1 = both_sources[['upid', 'biosample', 'protein_count_up', 'protein_count_atb']]

###############################################################################
# Preparing plot data for duplicate biosamples
df_up_biosample_dups = df_up_biosample_dups.sort_values(by=['biosample', 'protein_count'], ascending=[True, False])
df_up_bs_dups=df_up_biosample_dups[['upid', 'biosample', 'protein_count', 'complete_combined_score', 'is_excluded', 'is_effective']]
# Sample data
#data = {
#    'upid': ['UP000047540', 'UP000235454', 'UP000235499', 'UP000046519'],
#    'biosample': ['SAMEA1024557', 'SAMEA1024557', 'SAMEA1024588', 'SAMEA1024588'],
#    'protein_count': [2010, 1982, 2028, 2012],
#    'complete_combined_score': [99.8, 99.8, 99.8, 99.3],
#    'is_excluded': ['t', 'f', 'f', 'f'],
#    'is_effective': [True, False, False, False]
#}
#df = pd.DataFrame(data)

# Add a column to indicate row order within each biosample
df_up_bs_dups['pair_index'] = df_up_bs_dups.groupby('biosample').cumcount() + 1

# Reshape the data using pivot
df_pivot = df_up_bs_dups.pivot(index='biosample', columns='pair_index')

# Flatten the MultiIndex columns
df_pivot.columns = [f"{col}{num}" for col, num in df_pivot.columns]

# Reset index to make 'biosample' a column again
df_pivot = df_pivot.reset_index()

print(df_pivot)
print(f"\nChanged shape of data from: {df_up_bs_dups.shape} ---> to: {df_pivot.shape} ")

##################################
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Sample data
data = {
    "biosample": ["SAMEA1024557", "SAMEA1024588"],
    "upid1": ["UP000047540", "UP000235499"],
    "upid2": ["UP000235454", "UP000046519"],
    "protein_count1": [2010, 2028],
    "protein_count2": [1982, 2012],
    "complete_combined_score1": [99.8, 99.8],
    "complete_combined_score2": [99.8, 99.3],
}

# Create DataFrame
df = pd.DataFrame(data)

# Define positions and bar width
x = np.arange(len(df))  # Positions for biosamples on X-axis
bar_width = 0.2
gap = 0.3  # Gap between '1' and '2' groups

# Create the figure and axis objects
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot protein counts (ending in '1') on the left y-axis
ax1.bar(x - bar_width - gap / 2, df["protein_count1"], bar_width, label="Protein Count 1", color="blue")
ax1.bar(x - gap / 2, df["protein_count2"], bar_width, label="Protein Count 2", color="lightblue")
ax1.set_ylabel("Protein Count", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")

# Add text labels for protein counts
for i in range(len(df)):
    ax1.text(x[i] - bar_width - gap / 2, df["protein_count1"][i] + 10, str(df["protein_count1"][i]), ha="center", color="blue")
    ax1.text(x[i] - gap / 2, df["protein_count2"][i] + 10, str(df["protein_count2"][i]), ha="center", color="blue")

# Create a second y-axis for BUSCO scores (ending in '1' and '2')
ax2 = ax1.twinx()
ax2.bar(x + gap / 2, df["complete_combined_score1"], bar_width, label="BUSCO Score 1", color="black")
ax2.bar(x + bar_width + gap / 2, df["complete_combined_score2"], bar_width, label="BUSCO Score 2", color="gray")
ax2.set_ylabel("BUSCO Score", color="black")
ax2.tick_params(axis="y", labelcolor="black")

# Add text labels for BUSCO scores
for i in range(len(df)):
    ax2.text(x[i] + gap / 2, df["complete_combined_score1"][i] + 0.5, str(df["complete_combined_score1"][i]), ha="center", color="black")
    ax2.text(x[i] + bar_width + gap / 2, df["complete_combined_score2"][i] + 0.5, str(df["complete_combined_score2"][i]), ha="center", color="black")

# Set x-axis labels and title
plt.xticks(x, df["biosample"])
plt.title("Protein Counts and BUSCO Scores by Biosample")

# Add legend
fig.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=4)

# Show the plot
plt.tight_layout()
plt.show()
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Sample data
data = {
    "biosample": ["SAMEA1024557", "SAMEA1024588"],
    "upid1": ["UP000047540", "UP000235499"],
    "upid2": ["UP000235454", "UP000046519"],
    "protein_count1": [2010, 2028],
    "protein_count2": [1982, 2012],
    "complete_combined_score1": [99.8, 99.8],
    "complete_combined_score2": [99.8, 99.3],
    "is_excluded1": ["t", "f"],
    "is_excluded2": ["f", "f"],
    "is_effective1": [True, False],
    "is_effective2": [False, False],
}

# Create DataFrame
df = pd.DataFrame(data)

# Define positions and bar width
x = np.arange(len(df))  # Positions for biosamples on X-axis
bar_width = 0.15
gap = 0.15  # Slight gap between '1' and '2' groups

# Function to determine bar color based on conditions
def get_bar_color(is_excluded, is_effective):
    if is_excluded == "f" and not is_effective:
        return "green"
    elif is_excluded == "t" and is_effective:
        return "red"
    else:
        return "gold"

# Create the figure and axis objects
fig, ax1 = plt.subplots(figsize=(10, 6))

ax2 = ax1.twinx()

# Plot protein_count1 and complete_combined_score1 for each biosample (left group)
ax1.bar(
    x - bar_width - gap / 2,
    df["protein_count1"],
    bar_width,
    edgecolor = 'blue',
    label="Protein Count 1",
    color=[get_bar_color(df["is_excluded1"][i], df["is_effective1"][i]) for i in range(len(df))],
)
ax2.bar(
    x - gap / 2,
    df["complete_combined_score1"],
    bar_width,
    edgecolor = 'black',
    label="BUSCO Score 1",
    color=[get_bar_color(df["is_excluded1"][i], df["is_effective1"][i]) for i in range(len(df))],
)

# Plot protein_count2 and complete_combined_score2 for each biosample (right group)
ax1.bar(
    x + gap / 2,
    df["protein_count2"],
    bar_width,
    label="Protein Count 2",
    edgecolor = 'blue',
    color=[get_bar_color(df["is_excluded2"][i], df["is_effective2"][i]) for i in range(len(df))],
    hatch="//",  # Add hatches to differentiate columns ending in '2'
)
ax2.bar(
    x + bar_width + gap / 2,
    df["complete_combined_score2"],
    bar_width,
    label="BUSCO Score 2",
    edgecolor = 'black',
    color=[get_bar_color(df["is_excluded2"][i], df["is_effective2"][i]) for i in range(len(df))],
    hatch="//",  # Add hatches to differentiate columns ending in '2'
)

# Set labels for the left y-axis
ax1.set_ylabel("Protein Count", color="blue")
ax2.set_ylabel("BUSCO Score", color="black")

ax1.tick_params(axis="y", labelcolor="blue")
ax2.tick_params(axis="y", labelcolor="blue")


# Add text labels above bars for protein counts and BUSCO scores
for i in range(len(df)):
    # Protein Count 1
    ax1.text(x[i] - bar_width - gap / 2, df["protein_count1"][i] + 10, str(df["protein_count1"][i]), ha="center", color="blue")
    
    # BUSCO Score 1
    ax2.text(x[i] - gap / 2, df["complete_combined_score1"][i] + 0.5, str(df["complete_combined_score1"][i]), ha="center", color="black")
    
    # Protein Count 2
    ax1.text(x[i] + gap / 2, df["protein_count2"][i] + 10, str(df["protein_count2"][i]), ha="center", color="blue")
    
    # BUSCO Score 2
    ax2.text(x[i] + bar_width + gap / 2, df["complete_combined_score2"][i] + 0.5, str(df["complete_combined_score2"][i]), ha="center", color="black")

# Set x-axis labels and title
plt.xticks(x, df["biosample"])
plt.title("Protein Counts and BUSCO Scores by Biosample")

# Add legend
fig.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=4)

# Show the plot
plt.tight_layout()
plt.show()


###############################################################################
#===========
# Bar plot
#===========
# --- PLOT: matplotlib ---
categories = ["Only UP", "Both", "Only ATB"]
counts = [len(only_up), len(both_sources), len(only_atb)]

plt.figure(figsize=(8, 7))
plt.bar(categories, counts, color=["blue", "brown", "darkorange"], alpha=1)

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Proteomes")
plt.title("Comparison of Proteome counts by Biosample in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)both_sources

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
bar_plot = sns.barplot(x='Data Sources', y='Number of Proteomes', data=data, alpha=0.99, palette=["blue", "brown", "darkorange"])

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
venn.get_patch_by_id("01").set_color("darkorange")
venn.get_patch_by_id("11").set_color("brown")

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
            color=["blue", "brown", "darkorange"],
            alpha=1, 
            #color=["lightgrey", "lightgrey", "lightgrey"],
            width=bar_width)
# Overlay bars (colored with denser hatching for exclusive)
axes[1].bar(categories, 
            subset_counts, 
            color=["blue", "brown", "darkorange"], 
            alpha=1, 
            width=bar_width, 
            hatch='/')  # Denser hatching

# Labels & Title
axes[1].set_xlabel("Data Sources")
axes[1].set_ylabel("Number of Proteomes")
#axes[1].set_title("Stacked Bar Plot of Proteome Counts")
axes[1].grid(axis="y", 
             linestyle="--", 
             alpha=1)

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
    Patch(facecolor='lightgrey', hatch='/', label='Exclusive')  # Denser hatch
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
        color="orange", 
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
             alpha=0.8, 
             #element="step", 
             edgecolor='black', 
             linewidth=1)

# Histogram for ATB data
sns.histplot(data=df_atb2, 
             x='protein_count', 
             bins=bins, 
             color='darkorange',
             #label=f'ATB (n={len(df_atb2)}, Bin size: {bin_width}, Bin count: {num_bins_atb})',
             label=f'ATB (n={len(df_atb2)}, Bin count: {num_bins_atb})',
             alpha=0.7,
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
