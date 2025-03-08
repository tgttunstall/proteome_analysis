#!bin/env python
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import seaborn as sns
import itertool

#FIXME
#from plot_functions.py import plot_category_counts()
# also  needs data from plotcomp_atb_up_v1.py
###############################################################################
#===============
# Venn Diagram
#===============
# Define set sizes
only_up_count   = len(only_up)  # Unique to UP
only_atb_count  = len(only_atb)  # Unique to ATB
both_count      = len(both_sources)  # In both UP and ATB

# --- Venn Diagram ---
venn2(subsets= (only_up_count, only_atb_count, both_count), 
      set_labels=("UP",
                  "ATB"), 
      #ax=axes[0],
      set_colors=("blue", 
                  "darkorange"),
      alpha=0.7)
plt.title("Comparing Biosamples in UP and ATB")
plt.show()

#===============
# Stacked Bar plot Diagram
#===============
# Biosamples df
tot_up = only_up_count + both_count
tot_atb = only_atb_count + both_count
tot_up_atb = only_up_count + only_atb_count + both_count

data_biosamples = {'Categories': ['UP', 'ATB', 'UP+ATB'],
        'Total': [tot_up, tot_atb, tot_up_atb],
        'Exclusive': [only_up_count, only_atb_count, both_count]}

count_biosamples_df = pd.DataFrame(data_biosamples)
print(count_biosamples_df)

# Proteomes DF
dup_biosamples = 97
tot_up_proteomes     = tot_up + dup_biosamples
tot_atb_proteomes    = tot_atb
tot_up_atb_proteomes = tot_up_proteomes + tot_atb_proteomes
data_proteomes = {'Categories': ['UP', 'ATB', 'UP+ATB'],
        'Total': [tot_up_proteomes, tot_atb_proteomes, tot_up_atb_proteomes],
        'Exclusive': [only_up_count, only_atb_count, both_count]}

count_proteomes_df = pd.DataFrame(data_proteomes)
print(count_proteomes_df)

#--------
# Barplot: Biosamples
#--------
outplot_biosamples = basedir +"/Plots/UP_ATB_BiosampleCount.png"
plot_category_counts(df = count_biosamples_df,
                         category_col="Categories", 
                         total_col="Total", 
                         total_col_label = "Total",
                         exclusive_col="Exclusive",
                         exclusive_col_label = "Exclusive",
                         cat_colours = {"UP": "blue",
                                   "ATB": "darkorange",
                                   "UP+ATB": "purple"},
                         bar_annot_size = 12,
                         xlabel = "",
                         ylabel = "Count",
                         label_size = 14,
                         plot_title="Biosamples",
                         plot_label_size=16,
                         legend_font_size=12,
                         output_plot=outplot_biosamples)
                        

#--------
# Barplot: Proteomes
#--------
outplot_proteomes = basedir +"/Plots/UP_ATB_ProteomesCount.png"
plot_category_counts(df = count_proteomes_df,
                         category_col="Categories", 
                         total_col="Total", 
                         total_col_label = "Total",
                         exclusive_col="Exclusive",
                         exclusive_col_label = "Exclusive",
                         cat_colours = {"UP": "blue",
                                   "ATB": "darkorange",
                                   "UP+ATB": "purple"},
                         bar_annot_size = 12,
                         xlabel = "",
                         ylabel = "Count",
                         label_size = 14,
                         plot_title="Proteomes",
                         plot_label_size=16,
                         legend_font_size=12,
                         output_plot=outplot_proteomes)

###############################################################################
a = df_up2[['upid', 'biosample', 'is_excluded', 'is_effective','exclusion_id']]
not_true_exclusions = [96,99]
a.loc[a['exclusion_id'].str.contains(r'96|99'), :].shape[0]
a2 = a.loc[a['exclusion_id'].str.contains(r'96|99'), :]


nbins = np.histogram_bin_edges(cluster_sizes, bins=hist_bin_method)

# Combining 'Representative' and 'Reference' data into one Series
combined_data_busco = pd.concat([
    df_up2[df_up2['is_representative'] == 't']['complete_combined_score'], 
    df_up2[df_up2['is_reference'] == 't']['complete_combined_score']
])

red_df = df_up2[df_up2['is_redundant'].isin([1, -1])]
ex_df = df_up2[df_up2['is_excluded'] == 't']


df_up2[df_up2['is_excluded'] == 't']
df_up2[df_up2['is_redundant'].isin([1, -1])]

ax=sns.histplot(
    data=ex_df,
    x='checkm-completeness', 
    bins=50, 
    color='grey',
    label='', 
    alpha=excluded_alpha)
plt.yscale('log')

df_up2[df_up2['is_redundant'].isin([1, -1])]['complete_combined_score'].describe()

