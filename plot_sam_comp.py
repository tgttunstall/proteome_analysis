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
