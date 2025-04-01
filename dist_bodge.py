
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


data_df = df_up_all.copy()
bin_method='rice'
nbins = np.histogram_bin_edges(data_df['protein_count'], bins=bin_method)
num_bins = len(nbins) - 1
bin_size = round(nbins[1] - nbins[0], 2) if num_bins > 1 else "Variable"

#bins = np.arange(data_df['protein_count'].min(), data_df['protein_count'].max() + 100, 100)

title_fontsize = 11
label_fontsize = 14
tick_fontsize = 12
legend_fontsize = 10

general_alpha = 0.9
reference_alpha = 0.8
rep_alpha = 0.8
redundant_alpha = 0.5
excluded_alpha = 0.6
atb_alpha = 0.8

final_count = len(data_df)

# Create a 2x3 grid of subplots
fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharex=True, sharey=True)
axes = axes.flatten()  # Flatten to simplify indexing

# Plot configurations for each category
# General UP Proteomes with Reference and Representative overlays
sns.histplot(data=data_df, x='protein_count', 
             bins=num_bins, 
             color='blue', label='General UP Proteomes',
             alpha=general_alpha,
             ax=axes[0]
             )
sns.histplot(data=data_df[data_df['is_reference'] == 't'], x='protein_count', 
             bins=num_bins, 
             color='magenta',
             label='Reference', 
             alpha=reference_alpha, 
             hatch="--",
             ax=axes[0]
             )
sns.histplot(data=data_df[data_df['is_representative'] == 't'], x='protein_count', 
             bins=num_bins, 
             color='cyan',
             label='Representative', 
             alpha=rep_alpha, 
             hatch="////",
             ax=axes[0]
             )

#plt.yscale('log')
#plt.title(f'Distribution of All UP proteomes(n={final_count:,}) | Bins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
#plt.legend(fontsize=legend_fontsize)
#plt.show()


axes[0].set_yscale('log')
axes[0].set_title(f'Distribution of All UP proteomes(n={final_count:,}) | Bins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
#axes[0].legend(fontsize=legend_fontsize)
plt.show()


# New Detailed UP Categories Plot
sns.histplot(data=data_df, x='protein_count', bins=num_bins, 
             color='blue', 
             label='General UP Proteomes', 
             alpha=0.2, 
             ax=axes[0])
sns.histplot(data=data_df[data_df['is_redundant'].isin([1, -1])], x='protein_count', 
             bins=num_bins, 
             color='black', 
             label='Redundant', 
             alpha=0.3, 
             ax=axes[0])
sns.histplot(data=data_df[data_df['is_excluded'] == 't'], x='protein_count', 
             bins=num_bins, 
             color='red',
             label='Excluded', 
             alpha=0.3, 
             ax=axes[0])
#plt.yscale('log')
#plt.title('Detailed UP Categories', fontsize=title_fontsize)
#plt.legend(fontsize=legend_fontsize)
#plt.tight_layout()

axes[1].set_yscale('log')
axes[1].set_title(f'Distribution of All UP proteomes(n={final_count:,}) | Bins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
#axes[1].legend(fontsize=legend_fontsize)
plt.tight_layout()
plt.show()

##################################

# New Detailed UP Categories Plot
sns.histplot(data=df_up2, x='protein_count', bins=bins, 
             color='blue', 
             label='General UP Proteomes', 
             alpha=0.2, 
             #ax=axes[5]
             )
sns.histplot(data=df_up2[df_up2['is_redundant'].isin([1, -1])], x='protein_count', bins=bins, 
             color='black', 
             label='Redundant', 
             alpha=0.3, 
             #ax=axes[5]
             )
sns.histplot(data=df_up2[df_up2['is_excluded'] == 't'], x='protein_count', bins=bins, 
             color='red',
             label='Excluded', 
             alpha=0.3, 
             #ax=axes[5]
             )
plt.yscale('log')
plt.title(f'Distribution of All UP proteomes by category (n={final_count:,}) | Bins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
plt.legend(fontsize=legend_fontsize)

plt.tight_layout()
plt.show()

data_df['is_excluded'].value_counts()
data_df['is_redundant'].value_counts()
#data_df(['is_excluded', 'is_redundant']).value_counts()
pc = list(df_up_all['protein_count'])
pcC = pc.count(0)
pcC
df_up_all['protein_count'].nunique()

dfPC = df_up_all[df_up_all['protein_count']== 0]
dfPC['upid'].nunique()


###############################################################################
# UP and ATB
bin_method='rice'
nbins_up = np.histogram_bin_edges(df_up2['protein_count'], bins=bin_method)
num_bins = len(nbins) - 1
bin_size = round(nbins[1] - nbins[0], 2) if num_bins > 1 else "Variable"

# Calculate mean and standard deviation for both datasets
mean_up = df_up['protein_count'].mean()
std_up = df_up['protein_count'].std()

mean_atb = df_atb2['protein_count'].mean()
std_atb = df_atb2['protein_count'].std()

sns.histplot(data=df_up, 
             x='protein_count', 
             #bins=num_bins,
             binwidth=50,
             color='blue', 
             label='UP Proteomes', 
             alpha=0.8, 
             #ax=axes[5]
             )
sns.histplot(data=df_atb2,
             x='protein_count', 
             #bins=num_bins,
             binwidth=50,
             color='darkorange', 
             label='ATB Proteomes', 
             alpha=0.8, 
             #ax=axes[5]
             )
plt.yscale('log')
plt.xlabel("Number of Proteins")
plt.ylabel("Proteome Frequency (Log)")
plt.title(f'Distribution of UP and ATB proteomes | Bins: {num_bins}, Bin Size: {bin_size}', fontsize=title_fontsize)
plt.legend(fontsize=legend_fontsize)


# Add vertical lines for means
plt.axvline(mean_up, color='blue', linestyle='dashed', linewidth=2, label=f'UP Mean: {mean_up:.1f}')
plt.axvline(mean_atb, color='darkorange', linestyle='dashed', linewidth=2, label=f'ATB Mean: {mean_atb:.1f}')

# Add shading for ±1 SD
plt.fill_betweenx(y=[0, plt.gca().get_ylim()[1]], 
                  x1=mean_up - std_up, x2=mean_up + std_up, 
                  color='blue', alpha=0.2, label=f'UP ±1 SD')

plt.fill_betweenx(y=[0, plt.gca().get_ylim()[1]], 
                  x1=mean_atb - std_atb, x2=mean_atb + std_atb, 
                  color='darkorange', alpha=0.2, label=f'ATB ±1 SD')

# Annotate mean and SD values with manually adjusted positions
# plt.text(mean_up + 200, plt.gca().get_ylim()[1] * 0.3, 
#          f"Mean: {mean_up:.1f}\nSD: {std_up:.1f}", 
#          color='blue', fontsize=10,
#          #bbox=dict(facecolor='white', alpha=0.6)
#          )

# plt.text(mean_atb - 100, plt.gca().get_ylim()[1] * 0.5,
#          color='darkorange', fontsize=10, 
#          #bbox=dict(facecolor='white', alpha=0.6)
#          )


# Add mean and SD info just below the legend
plt.figtext(0.15, 0.6, f"Mean = {mean_up:.1f}, \nSD = {std_up:.1f}", 
            fontsize=10, color='blue', ha='left')
plt.figtext(0.15, 0.52, f"Mean = {mean_atb:.1f}, \nSD = {std_atb:.1f}", 
            fontsize=10, color='darkorange', ha='left')

plt.tight_layout()
plt.show()
###############################################################################

data_df['is_excluded'].value_counts()
data_df['is_redundant'].value_counts()
#data_df(['is_excluded', 'is_redundant']).value_counts()
pc = list(df_up_all['protein_count'])
pcC = pc.count(0)
pcC
df_up_all['protein_count'].nunique()

dfPC = df_up_all[df_up_all['protein_count']== 0]
dfPC['upid'].nunique()

