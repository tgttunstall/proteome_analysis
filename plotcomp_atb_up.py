#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# =============================================================================
# # Generate dummy data
# def generate_sample_id(prefix, num):
#     return f"{prefix}{str(num).zfill(6)}"
# 
# # Generate all samples
# num_all_samples = 10000
# all_samples = pd.DataFrame({
#     'sample_id': [generate_sample_id('SAM', i) for i in range(num_all_samples)],
#     'protein_count': np.random.normal(2500, 300, num_all_samples).astype(int)
# })
# 
# # Generate mapped samples (subset of all samples)
# num_mapped_samples = 7000
# mapped_samples = all_samples.sample(n=num_mapped_samples).copy()
# mapped_samples['upid'] = [f"UP{str(i).zfill(6)}" for i in range(num_mapped_samples)]
# 
# # Save generated data to TSV files
# all_samples.to_csv('all_samples.tsv', sep='\t', index=False, header=False)
# mapped_samples.to_csv('mapped_samples.tsv', sep='\t', index=False, header=False)
# =============================================================================
#####
# R2
#####

homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

# Load TSV files
#up_file = basedir + "/up_pcounts.tsv"
up_file = basedir + "/spneumo_biosample_info.out"
atb_file = basedir + "/atb_pcounts.tsv"

# Read the data
#df_up = pd.read_csv(up_file, sep="\t")
df_up = pd.read_csv(up_file, sep="\t", usecols = ['biosample', 'protein_count'])
df_atb = pd.read_csv(atb_file, sep="\t")

# sanity check: Drop duplicated, and data with protein count == 0

df_up2 = df_up[df_up['protein_count'] > 0]
df_atb2 = df_atb[df_atb['protein_count'] > 0]

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
# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(df_merged["protein_count_up"], bins=50, alpha=0.5, label="UP")
plt.hist(df_merged["protein_count_atb"], bins=50, alpha=0.5, label="ATB")
plt.xlabel("Protein Count")
plt.ylabel("Frequency")
plt.title("Protein Count Distribution Across Data Sources")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()

###############################################################################
#===========
# Bar plot
#===========
categories = ["Only UP", "Only ATB", "Both"]
counts = [len(only_up), len(only_atb), len(both_sources)]

# Plot
plt.figure(figsize=(8, 7))
plt.bar(categories, counts, color=["blue", "red", "green"], alpha=0.7)

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Proteomes")
plt.title("Comparison of Proteome counts by Biosample in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts
for i, v in enumerate(counts):
    plt.text(i, v + max(counts) * 0.02, str(v), ha="center", fontsize=14)

plt.show()

###############################################################################
#================================================
# Stacked bar plot with total for each category
#================================================
categories = ["UP", "UP+ATB", "ATB"]
total_counts = [len(only_up) + len(both_sources),  # Total UP
                len(only_up) + len(only_atb) + len(both_sources),# Total UP+ATB
                len(only_atb) + len(both_sources)] # Total ATB]  
subset_counts = [len(only_up), len(both_sources), len(only_atb)]  # Stacked portion

# Plot
plt.figure(figsize=(8, 7))
bar_width = 0.6

# Base bars (lighter colors, transparent)
#plt.bar(categories, total_counts, color=["lightblue", "lightgreen", "lightcoral"], alpha=0.5, label="Total", width=bar_width)
plt.bar(categories, total_counts, color=["lightgrey", "lightgrey", "lightgrey"], 
        alpha=0.5
        , label="Total", 
        width=bar_width)

# Overlay bars (darker, non-transparent for actual subset counts)
plt.bar(categories, subset_counts, color=["blue", "green", "red"], 
        alpha=1, 
        label="Subset", 
        width=bar_width)

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Proteomes")
plt.title("Stacked Bar Plot of Proteome Counts in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate total and subset counts
for i, (total, subset) in enumerate(zip(total_counts, subset_counts)):
    plt.text(i, 
             total + max(total_counts) * 0.008, 
             str(total), 
             ha="center", 
             fontsize=12, 
             color="black")
    plt.text(i, 
             subset / 2, 
             str(subset), 
             ha="center", 
             fontsize=12, 
             color="white", 
             fontweight="bold")  # Centered inside subset bar

plt.legend(["Total", "Exclusive"])
plt.show()
###############################################################################
#======
# Venn
#======
from matplotlib_venn import venn2

# Define set sizes
only_up_count = len(only_up)  # Unique to UP
only_atb_count = len(only_atb)  # Unique to ATB
both_count = len(both_sources)  # In both UP and ATB

# Create Venn diagram
plt.figure(figsize=(6, 6))
venn = venn2(subsets=(only_up_count, 
                      only_atb_count, 
                      both_count), set_labels=("UP", "ATB"))

# Customize text labels
venn.get_label_by_id("10").set_text(f"{only_up_count}")  # Only UP
venn.get_label_by_id("01").set_text(f"{only_atb_count}")  # Only ATB
venn.get_label_by_id("11").set_text(f"{both_count}")  # Both UP & ATB

# Add title
plt.title("Comparing Proteomes/Biosamples in UP and ATB")
plt.show()

###############################################################################
#===============
# side by side: venn and BP
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

# Base bars (lighter colors)
axes[1].bar(categories, 
            total_counts, 
            #color=["lightblue", "lightcoral", "lightgreen"], 
            color=["lightgrey", "lightgrey", "lightgrey"], 
            alpha=0.5, 
            width=bar_width, 
            label="Total")

# Overlay bars (darker colors for actual subset counts)
axes[1].bar(categories, 
            subset_counts, 
            color=["blue", "green", "red"], 
            alpha=1, 
            width=bar_width, 
            label="Subset")

# Labels & Title
axes[1].set_xlabel("Data Sources")
axes[1].set_ylabel("Number of Proteomes")
axes[1].set_title("Stacked Bar Plot of Proteome Counts")
axes[1].grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts on bars
for i, (total, subset) in enumerate(zip(total_counts, subset_counts)):
    axes[1].text(i, 
                 total + max(total_counts) * 0.008, 
                 str(total), 
                 ha="center", 
                 fontsize=12, 
                 color="black")
    axes[1].text(i, subset / 2, str(subset), 
                 ha="center", 
                 fontsize=12, 
                 color="white", 
                 fontweight="bold")  # Centered inside subset bar

axes[1].legend(["Total", "Exclusive"])

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

# Plot histogram
plt.figure(figsize=(8, 6))
plt.hist(up_proteins, bins=50, alpha=0.6, color="blue", label="UP", edgecolor="black")
plt.hist(atb_proteins, bins=50, alpha=0.6, color="red", label="ATB", edgecolor="black")

plt.yscale("log")

# Labels & Title
plt.xlabel("Number of Proteins")
plt.ylabel("Number of Proteomes")
plt.title("Protein Count Distribution Across Proteomes")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Show plot
plt.show()

###############################################################################
# Define bins (e.g., group by protein count ranges)
bins = np.arange(1, max(up_proteins.max(), atb_proteins.max()) + 100, 100)  # Adjust bin size as needed

# Count occurrences per bin
up_counts, _ = np.histogram(up_proteins, bins)
atb_counts, _ = np.histogram(atb_proteins, bins)

# Bar plot
plt.figure(figsize=(8, 6))
plt.bar(bins[:-1], up_counts, width=100, alpha=0.6, color="blue", label="UP", edgecolor="black")
plt.bar(bins[:-1], atb_counts, width=100, alpha=0.6, color="red", label="ATB", edgecolor="black")

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

###############################################################################
# protein numbers
# Assuming df_up and df_atb are your dataframes and 'protein_count' is the column of interest
total_proteins_up = up_proteins.sum()
average_proteins_up = up_proteins.mean()

total_proteins_atb = atb_proteins.sum()
average_proteins_atb = atb_proteins.mean()

data = {
    'Total Proteins': [total_proteins_up, total_proteins_atb],
    'Average Proteins': [average_proteins_up, average_proteins_atb]
}
stats_df = pd.DataFrame(data, index=['UP', 'ATB'])


import matplotlib.pyplot as plt

# Plotting
fig, ax1 = plt.subplots(figsize=(10, 6))

# Bar plot for total proteins
totals = ax1.bar(stats_df.index, stats_df['Total Proteins'], color='b', alpha=0.6, label='Total Proteins')
# Labeling total proteins
for rect in totals:
    height = rect.get_height()
    ax1.annotate('{}'.format(int(height)),
                 xy=(rect.get_x() + rect.get_width() / 2, height),
                 xytext=(0, 3),  # 3 points vertical offset
                 textcoords="offset points",
                 ha='center', va='bottom')

# Create a second y-axis for the averages
ax2 = ax1.twinx()
averages = ax2.bar(stats_df.index, stats_df['Average Proteins'], color='r', alpha=0.6, label='Average Proteins', width=0.4)
# Labeling average proteins
for rect in averages:
    height = rect.get_height()
    ax2.annotate('{:.2f}'.format(height),
                 xy=(rect.get_x() + rect.get_width() / 2, height),
                 xytext=(0, 3),  # 3 points vertical offset
                 textcoords="offset points",
                 ha='center', va='bottom')

# Setting labels and titles
ax1.set_xlabel('Dataset', fontsize=12)
ax1.set_ylabel('Total Protein Count', color='b', fontsize=12)
ax2.set_ylabel('Average Protein Count', color='r', fontsize=12)
ax1.set_title('Comparison of Total and Average Protein Counts in UP and ATB Datasets', fontsize=14)

# Adding legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper left')

# Show the plot
plt.show()


##############################################################################
import pandas as pd

# Load the data from a CSV file
df_up_all = pd.read_csv(up_file, delimiter='\t')
df_up_all2 = df_up_all[df_up_all['protein_count'] > 0]
df_up_all2 = df_up_all.drop_duplicates(subset = 'biosample', keep = 'last')


# count values
# Count the values for each relevant column
representative_count = df_up_all2['is_representative'].value_counts()
reference_count = df_up_all2['is_reference'].value_counts()
redundant_count = df_up_all2['is_redundant'].value_counts()
excluded_count = df_up_all2['is_excluded'].value_counts()

# Prepare a DataFrame from the counts
status_data = pd.DataFrame({
    'Representative': representative_count,
    'Reference': reference_count,
    #'Redundant': redundant_count,
    'Excluded': excluded_count
}).fillna(0)  # Fill NaNs with 0 in case some categories don't exist in the data

################################################################################
import matplotlib.pyplot as plt
import numpy as np

# Extract counts
up_proteins = df_up_all2["protein_count"]
atb_proteins = df_up_all2["protein_count"]

# Create masks for each category
is_excluded = (df_up_all2["is_excluded"] == 1) | (df_up2["is_excluded"] == -1)
is_redundant = df_up_all2["is_redundant"] == 't'
is_reference = df_up_all2["is_reference"] == 't'
is_representative = df_up_all2["is_representative"] == 't'

# Plot histogram
plt.figure(figsize=(12, 8))

# Plot ATB data
plt.hist(atb_proteins, bins=50, alpha=0.6, color="red", label="ATB", edgecolor="black")

# Plot UP data with different categories, excluding those already plotted
up_unclassified = ~(is_excluded | is_redundant | is_reference | is_representative)
plt.hist(up_proteins[up_unclassified], bins=50, alpha=0.6, color="blue", 
         label="UP (Other)", edgecolor="black")

up_excluded = is_excluded & up_unclassified
plt.hist(up_proteins[up_excluded], bins=50, alpha=0.6, color="yellow", 
         label="UP (Excluded)", edgecolor="black")

up_redundant = is_redundant & up_unclassified
plt.hist(up_proteins[up_redundant], bins=50, alpha=0.6, color="orange", 
         label="UP (Redundant)", edgecolor="black")

up_reference = is_reference & up_unclassified
plt.hist(up_proteins[up_reference], bins=50, alpha=0.6, color="purple", 
         label="UP (Reference)", edgecolor="black")

up_representative = is_representative & up_unclassified
plt.hist(up_proteins[up_representative], bins=50, alpha=0.6, color="green", 
         label="UP (Representative)", edgecolor="black")

plt.yscale("log")

# Labels & Title
plt.xlabel("Number of Proteins")
plt.ylabel("Number of Proteomes")
plt.title("Protein Count Distribution Across Proteomes")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Show plot
plt.show()


###############################################################################

#####
# R1
#####
# =============================================================================
# Read the data 
all_samples = pd.read_csv('/home/pub/Work/data_arise_proteome/spneumo_dataset/del/p_counts.txt', 
                          sep='\t', 
                          names=['sample_id', 'protein_count'])
mapped_samples = pd.read_csv('/home/pub/Work/data_arise_proteome/spneumo_dataset/del/SAM_proteins_per_UPID.tsv',
                             sep='\t',
                             skiprows=1,
                             names=['upid', 'sample_id', 'protein_count'])



# Create the plot
plt.figure(figsize=(12, 6))

# Plot histograms
plt.hist(all_samples['protein_count'], bins=50, alpha=0.5, label='All Samples')
plt.hist(mapped_samples['protein_count'], bins=50, alpha=0.5, label='Mapped Samples')

# Customize the plot
plt.title('Distribution of Protein Counts: All Samples vs Mapped Samples', fontsize=16)
plt.xlabel('Number of Proteins', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)

## v2
plt.figure(figsize=(12, 6))

# Plot histograms
plt.hist(all_samples['protein_count'], bins=40, range=(0, 4000), alpha=0.5, label='All Samples')
plt.hist(mapped_samples['protein_count'], bins=40, range=(0, 4000), alpha=0.5, label='Mapped Samples')

# Customize the plot
plt.title('Distribution of Protein Counts: All Samples vs Mapped Samples', fontsize=16)
plt.xlabel('Number of Proteins', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xlim(1000,3000)  # Set x-axis limits explicitly
plt.legend()
plt.grid(True, alpha=0.3)

# # Add statistics
# for df, label in [(all_samples, 'All'), (mapped_samples, 'Mapped')]:
#     mean = df['protein_count'].mean()
#     median = df['protein_count'].median()
#     plt.axvline(mean, color='r' if label == 'All' else 'g', linestyle='dashed', 
#                 linewidth=2, label=f'{label} Mean: {mean:.2f}')
#     plt.axvline(median, color='r' if label == 'All' else 'g', linestyle='dotted', 
#                 linewidth=2, label=f'{label} Median: {median:.2f}')

# Show the plot
plt.legend()
plt.tight_layout()
plt.show()

# Print additional statistics
print(f"Total samples: {len(all_samples)}")
print(f"Mapped samples: {len(mapped_samples)}")
print(f"Percentage mapped: {len(mapped_samples)/len(all_samples)*100:.2f}%")
###############################################################################
# Checking data
up_ds = homedir + "/Documents/arise/spneumo_dataset/spneumo_biosample_info.out"
atb_ds = homedir + "/Documents/arise/spneumo_dataset/atb_pcounts.tsv"

up_dataset = pd.read_csv(up_ds, sep = '\t')

up_dataset_corrupt = up_dataset[up_dataset['protein_count'] < 25]
up_dataset_corrupt = up_dataset[up_dataset['protein_count'] == 0]

corrupt_ds_name =  homedir + "/Documents/arise/spneumo_dataset/not_in_up.tsv"
up_dataset_corrupt.to_csv(corrupt_ds_name, sep = "\t")
