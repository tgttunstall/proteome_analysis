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
import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


homedir = os.path.expanduser("~")
# Load TSV files
#up_file = homedir + "/Documents/arise/spneumo_dataset/up_pcounts.tsv"
up_file = homedir + "/Documents/arise/spneumo_dataset/spneumo_biosample_info_clean.out"
atb_file = homedir + "/Documents/arise/spneumo_dataset/atb_pcounts.tsv"

# Read the data
#df_up = pd.read_csv(up_file, sep="\t")
df_up = pd.read_csv(up_file, sep="\t", usecols = ['biosample', 'protein_count'])
df_atb = pd.read_csv(atb_file, sep="\t")

# sanity check: Drop duplicated:
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
# diff type of plot
# Prepare stacked bar plot
categories = ["Only UP", "Only ATB", "Both"]
counts = [len(only_up), len(only_atb), len(both_sources)]

# Plot
plt.figure(figsize=(8, 7))
plt.bar(categories, counts, color=["blue", "red", "green"], alpha=0.7)

# Labels & Title
plt.xlabel("Data Sources")
plt.ylabel("Number of Biosamples")
plt.title("Comparison of Proteome counts by Biosample in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts
for i, v in enumerate(counts):
    plt.text(i, v + max(counts) * 0.02, str(v), ha="center", fontsize=14)

plt.show()

###############
# Data
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
plt.ylabel("Number of Biosamples")
plt.title("Stacked Bar Plot of Proteome Counts in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate total and subset counts
for i, (total, subset) in enumerate(zip(total_counts, subset_counts)):
    plt.text(i, 
             total + max(total_counts) * 0.02, 
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

plt.legend(["Total", "Subset"])
plt.show()

# Venn
from matplotlib_venn import venn2

# Define set sizes
only_up_count = len(only_up)  # Unique to UP
only_atb_count = len(only_atb)  # Unique to ATB
both_count = len(both_sources)  # In both UP and ATB

# Create Venn diagram
plt.figure(figsize=(6, 6))
venn = venn2(subsets=(only_up_count, only_atb_count, both_count), set_labels=("UP", "ATB"))

# Customize text labels
venn.get_label_by_id("10").set_text(f"{only_up_count}")  # Only UP
venn.get_label_by_id("01").set_text(f"{only_atb_count}")  # Only ATB
venn.get_label_by_id("11").set_text(f"{both_count}")  # Both UP & ATB

# Add title
plt.title("Overlap of Biosamples in UP and ATB")
plt.show()

###############################################################################
# side by side
# Define set sizes
only_up_count = len(only_up)  # Unique to UP
only_atb_count = len(only_atb)  # Unique to ATB
both_count = len(both_sources)  # In both UP and ATB

# Define bar plot data
categories = ["UP", "ATB", "UP+ATB"]
total_counts = [only_up_count + both_count, 
                only_up_count + only_atb_count + both_count,
                only_atb_count + both_count]
subset_counts = [only_up_count, only_atb_count, both_count]

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

axes[0].set_title("Overlap of Biosamples in UP and ATB")

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
            color=["blue", "red", "green"], 
            alpha=1, 
            width=bar_width, 
            label="Subset")

# Labels & Title
axes[1].set_xlabel("Data Sources")
axes[1].set_ylabel("Number of Biosamples")
axes[1].set_title("Stacked Bar Plot of Proteome Counts")
axes[1].grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts on bars
for i, (total, subset) in enumerate(zip(total_counts, subset_counts)):
    axes[1].text(i, 
                 total + max(total_counts) * 0.02, 
                 str(total), 
                 ha="center", 
                 fontsize=12, 
                 color="black")
    axes[1].text(i, subset / 2, str(subset), 
                 ha="center", 
                 fontsize=12, 
                 color="white", 
                 fontweight="bold")  # Centered inside subset bar

axes[1].legend(["Total", "Subset"])

# Show both plots
plt.tight_layout()
plt.show()

###############################################################################
# Extract protein counts
up_proteins = df_up["protein_count"]
atb_proteins = df_atb["protein_count"]

# Plot histogram
plt.figure(figsize=(8, 6))
plt.hist(up_proteins, bins=50, alpha=0.6, color="blue", label="UP", edgecolor="black")
plt.hist(atb_proteins, bins=50, alpha=0.6, color="red", label="ATB", edgecolor="black")

plt.yscale("log")

# Labels & Title
plt.xlabel("Number of Proteins per Biosample")
plt.ylabel("Number of Biosamples")
plt.title("Protein Count Distribution Across Biosamples")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Show plot
plt.show()

###############################################################################
import numpy as np

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
plt.xlabel("Number of Proteins per Biosample")
plt.ylabel("Number of Biosamples")
plt.title("Binned Protein Counts Across Biosamples")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)

plt.show()

###############################################################################

#####
# R1
#####
# =============================================================================
# Read the data 
all_samples = pd.read_csv('/home/pub/Work/data_arise_proteome/spneumo_dataset/p_counts.txt', sep='\t', names=['sample_id', 'protein_count'])
mapped_samples = pd.read_csv('/home/pub/Work/data_arise_proteome/spneumo_dataset/SAM_proteins_per_UPID.csv',
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
