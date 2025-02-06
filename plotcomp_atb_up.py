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
up_file = homedir + "/Documents/arise/spneumo_dataset/up_pcounts.tsv"
atb_file = homedir + "/Documents/arise/spneumo_dataset/atb_pcounts.tsv"

# Read the data
df_up = pd.read_csv(up_file, sep="\t")
df_atb = pd.read_csv(atb_file, sep="\t")

# Merge data on 'biosample' to see coverage
df_merged = pd.merge(df_up, df_atb, on="biosample", how="outer", suffixes=("_up", "_atb"))

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



# diff type of plot
# Prepare stacked bar plot
categories = ["Only UP", "Only ATB", "Both"]
counts = [len(only_up), len(only_atb), len(both_sources)]

# Plot
plt.figure(figsize=(8, 5))
plt.bar(categories, counts, color=["blue", "red", "green"], alpha=0.7)

# Labels & Title
plt.xlabel("Data Source Coverage")
plt.ylabel("Number of Biosamples")
plt.title("Comparison of Protein Coverage in UP and ATB")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Annotate counts
for i, v in enumerate(counts):
    plt.text(i, v + max(counts) * 0.02, str(v), ha="center", fontsize=12)

plt.show()


# diff plot version 2:
# Total dataset sizes
total_up = len(df_up)
total_atb = len(df_atb)

# Values for stacked bar chart
up_values = [len(only_up), len(both_sources)]  # UP breakdown (Only UP, Both)
#atb_values = [len(only_atb), len(both_sources)]  # ATB breakdown (Only ATB, Both)

up_values = [len(only_up), len(total_up)]  # UP breakdown (Only UP, Both)
atb_values = [len(only_atb), len(total_atb)]  # ATB breakdown (Only ATB, Both)


# Plot stacked bar chart
fig, ax = plt.subplots(figsize=(8, 6))

bar_width = 0.4  # Bar width for better spacing
x = np.arange(2)  # Two bars (UP and ATB)

# Stacked bars
ax.bar(x, [total_up, total_atb], color="lightgray", width=bar_width, label="Total")
ax.bar(x, up_values, color="blue", width=bar_width, label="Only UP or Shared")
ax.bar(x, atb_values, color="red", width=bar_width, bottom=up_values, label="Only ATB")

# Labels & Title
ax.set_xticks(x)
ax.set_xticklabels(["UP", "ATB"])
ax.set_ylabel("Number of Biosamples")
ax.set_title("Protein Coverage Breakdown in UP & ATB")
ax.legend()

# Annotate counts on bars
for i in range(2):  # Two bars
    ax.text(i, total_up if i == 0 else total_atb, str(total_up if i == 0 else total_atb), 
            ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.text(i, up_values[i], str(up_values[i]), ha="center", va="center", color="white", fontsize=11)
    ax.text(i, up_values[i] + atb_values[i], str(atb_values[i]), ha="center", va="center", color="white", fontsize=11)

plt.show()
   
###############################################################################
# Extract protein counts
up_proteins = df_up["protein_count"]
atb_proteins = df_atb["protein_count"]

# Plot histogram
plt.figure(figsize=(8, 6))
plt.hist(up_proteins, bins=50, alpha=0.6, color="blue", label="UP", edgecolor="black")
plt.hist(atb_proteins, bins=50, alpha=0.6, color="red", label="ATB", edgecolor="black")

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
bins = np.arange(0, max(up_proteins.max(), atb_proteins.max()) + 100, 100)  # Adjust bin size as needed

# Count occurrences per bin
up_counts, _ = np.histogram(up_proteins, bins)
atb_counts, _ = np.histogram(atb_proteins, bins)

# Bar plot
plt.figure(figsize=(8, 6))
plt.bar(bins[:-1], up_counts, width=100, alpha=0.6, color="blue", label="UP", edgecolor="black")
plt.bar(bins[:-1], atb_counts, width=100, alpha=0.6, color="red", label="ATB", edgecolor="black")

# Labels & Title
plt.xlabel("Number of Proteins per Biosample")
plt.ylabel("Number of Biosamples")
plt.title("Binned Protein Counts Across Biosamples")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.5)

plt.show()
###############################################################################
# Define x-axis limits
min_protein = min(up_proteins.min(), atb_proteins.min())  # Minimum protein count in both datasets
max_protein = max(up_proteins.max(), atb_proteins.max())  # Maximum protein count in both datasets
margin = 50  # Adjust margin as needed

x_min, x_max = min_protein - margin, max_protein + margin  # Expand range
bins = np.linspace(x_min, x_max, 50)  # 50 bins between min and max

# Plot histogram
plt.figure(figsize=(8, 6))
plt.hist(up_proteins, bins=bins, alpha=0.6, color="blue", label="UniProt (UP)", edgecolor="black")
plt.hist(atb_proteins, bins=bins, alpha=0.6, color="red", label="AntibioticDB (ATB)", edgecolor="black")

# Labels & Title
plt.xlabel("Number of Proteins per Biosample")
plt.ylabel("Number of Biosamples")
plt.title("Protein Count Distribution Across Biosamples")
plt.xlim(x_min, x_max)  # Apply x-axis range
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
