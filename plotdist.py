import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#sns.kdeplot(data=df_up2[df_up2['is_reference'] == 't']['protein_count'], color='magenta', label='Reference')
#sns.kdeplot(data=df_up2[df_up2['is_representative'] == 't']['protein_count'], color='cyan', label='Representative', linestyle="--", linewidth=10, warn_singular=False)

# Combining 'Representative' and 'Reference' data into one Series
combined_data = pd.concat([
    df_up2[df_up2['is_representative'] == 't']['protein_count'], 
    df_up2[df_up2['is_reference'] == 't']['protein_count']
])

plt.figure(figsize=(12, 8))

# Use seaborn.kdeplot to plot the combined data
#sns.kdeplot(data=combined_data, color='purple', label='Combined Rep & Ref', linewidth=2.5)

# Plot other categories
sns.kdeplot(data=df_up2['protein_count'], color='blue', label='General UP Proteomes', linewidth=2.5, linestyle="--")
sns.kdeplot(data=df_up2[df_up2['is_redundant'].isin([1, -1])]['protein_count'], color='black', label='Redundant', linewidth=2, alpha =0.5)
sns.kdeplot(data=df_up2[df_up2['is_excluded'] == 't']['protein_count'], color='red', label='Excluded')
sns.kdeplot(data=combined_data, color='magenta', label='Rep & Ref', linewidth=2.5, alpha = 0.5)
sns.kdeplot(data=df_atb2['protein_count'], color='darkorange', label='ATB Proteomes')


plt.xlabel('Number of Proteins', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.title('Proteome Distribution Across Categories', fontsize=16)
plt.legend(title='', fontsize=12)
plt.show()

####
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Assuming df_up2 and df_atb2 are preloaded DataFrames

# Bins and aesthetic configurations
bins = np.arange(df_up2['protein_count'].min(), df_up2['protein_count'].max() + 100, 100)
title_fontsize = 14
label_fontsize = 14
tick_fontsize = 12
legend_fontsize = 14

general_alpha = 0.9
reference_alpha = 0.8
rep_alpha = 0.8
redundant_alpha = 0.5
excluded_alpha = 0.6
atb_alpha = 0.8

# Create a 2x3 grid of subplots
fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
axes = axes.flatten()  # Flatten to simplify indexing

# Plot configurations for each category
# General UP Proteomes with Reference and Representative overlays
sns.histplot(data=df_up2, x='protein_count', bins=bins, 
             color='blue', label='General UP Proteomes',
             alpha=general_alpha,
             ax=axes[0])
sns.histplot(data=df_up2[df_up2['is_reference'] == 't'], x='protein_count', bins=bins, 
             color='magenta',
             label='Reference', 
             alpha=reference_alpha, 
             hatch="--",
             ax=axes[0])
sns.histplot(data=df_up2[df_up2['is_representative'] == 't'], x='protein_count', bins=bins, 
             color='cyan',
             label='Representative', 
             alpha=rep_alpha, 
             hatch="////",
             ax=axes[0])
axes[0].set_yscale('log')
axes[0].set_title('General UP Proteomes with Subcategories', fontsize=title_fontsize)
axes[0].legend(fontsize=legend_fontsize)

# Redundant
sns.histplot(data=df_up2[df_up2['is_redundant'].isin([1, -1])], x='protein_count', bins=bins, 
             color='black',
             label='Redundant', 
             alpha=redundant_alpha, 
             ax=axes[1])
axes[1].set_yscale('log')
axes[1].set_title('Redundant', fontsize=title_fontsize)
axes[1].legend(fontsize=legend_fontsize)

# Excluded
sns.histplot(data=df_up2[df_up2['is_excluded'] == 't'], x='protein_count', bins=bins, 
             color='red',
             label='Excluded', 
             alpha=excluded_alpha, 
             ax=axes[2])
axes[2].set_yscale('log')
axes[2].set_title('Excluded', fontsize=title_fontsize)
axes[2].legend(fontsize=legend_fontsize)

# ATB Proteomes
sns.histplot(data=df_atb2, x='protein_count', bins=bins, 
             color='darkorange', 
             label='ATB Proteomes', 
             alpha=atb_alpha, 
             ax=axes[3])
axes[3].set_yscale('log')
axes[3].set_title('ATB Proteome Count', fontsize=title_fontsize)
axes[3].legend(fontsize=legend_fontsize)

# All UP and ATB Proteomes
sns.histplot(data=df_up2, x='protein_count', bins=bins, 
             color='blue', label='All UP', 
             alpha=general_alpha-0.1, 
             ax=axes[4])
sns.histplot(data=df_atb2, x='protein_count', bins=bins, 
             color='darkorange', 
             label='All ATB', 
             alpha=atb_alpha-0.1, 
             ax=axes[4])
axes[4].set_yscale('log')
axes[4].set_title('All UP & ATB Proteomes', fontsize=title_fontsize)
axes[4].legend(fontsize=legend_fontsize)

# New Detailed UP Categories Plot
sns.histplot(data=df_up2, x='protein_count', bins=bins, 
             color='blue', 
             label='General UP Proteomes', 
             alpha=0.2, 
             ax=axes[5])
sns.histplot(data=df_up2[df_up2['is_redundant'].isin([1, -1])], x='protein_count', bins=bins, 
             color='black', 
             label='Redundant', 
             alpha=0.3, 
             ax=axes[5])
sns.histplot(data=df_up2[df_up2['is_excluded'] == 't'], x='protein_count', bins=bins, 
             color='red',
             label='Excluded', 
             alpha=0.3, 
             ax=axes[5])
axes[5].set_yscale('log')
axes[5].set_title('Detailed UP Categories', fontsize=title_fontsize)
axes[5].legend(fontsize=legend_fontsize)

plt.tight_layout()
plt.show()


##########
# GI
#########
#import plotly.graph_objects as go
#from cmplot import cmplot

#tips = sns.load_dataset("tips")
#iris = sns.load_dataset("iris")

#cmplot(tips,xcol='total_bill')
#fig = go.Figure(*cmplot(iris,xcol="species")) #using splat operator
#fig.show()
