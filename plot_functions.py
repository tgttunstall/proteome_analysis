#!/bin/env python

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

###############################################################################
def plot_clustersize(df, 
              plot_colname='proteome_count', 
              total_proteomes=None, 
              total_proteins=None,
              min_proteomes=None,
              hist_bin_method='rice',
              ylog=True,
              xlabel='Unique Proteomes',
              ylabel='Number of Clusters',
              x_proteomes='number', 
              output_plot='outplot.png', 
              show_stats=True):
    """
    Generates a histogram of cluster sizes from a precomputed counts file.

    Args:
        df: DataFrame containing the data.
        plot_colname (str): Column name to plot.
        total_proteomes (int): Total number of proteomes in the dataset (needed for % calculation).
        output_plot (str): Path to save the output plot file.
        min_proteomes (int): Minimum number of proteomes to filter by.
        ylog (bool): Whether to apply logarithmic scale on y-axis.
        x_proteomes (str): 'number' for raw counts, 'percent' for percentages.
        show_stats (bool): If True, prints and overlays core/accessory/rare stats on the plot.
    """
    # Check if the column name exists in the DataFrame
    if plot_colname not in df.columns:
        print(f"Column name '{plot_colname}' not found in DataFrame.")
        return
    
    print(f"\nLoading input data for plotting. Number of initial entries: {df.shape[0]}")

    # Apply filtering
    if min_proteomes is not None: 
        df_filtered = df[df[plot_colname] > min_proteomes].dropna(subset=[plot_colname])
        if df_filtered.empty:
            print("No clusters left after filtering based on minimum proteomes.")
            return
        print(f"Entries after filtering: {df_filtered.shape[0]}")
    else:
        df_filtered = df

    # Convert to percentages if required
    #FIXME: check if percent or number are there as only  these are valid params
    if x_proteomes == 'percent':
        if total_proteomes is None:
            print("Error: total_proteomes must be provided when x_proteomes='percent'")
            return
        cluster_sizes = (df_filtered[plot_colname] / total_proteomes) * 100
        nbins = np.histogram_bin_edges(cluster_sizes, bins=hist_bin_method)
        
    else:
        cluster_sizes = df_filtered[plot_colname]
        nbins = np.histogram_bin_edges(cluster_sizes, bins=hist_bin_method)

    print(f"\nNo. of bins with method {hist_bin_method}: {len(nbins)}")

    num_bins = len(nbins) - 1
    bin_size = round(nbins[1] - nbins[0], 2) if num_bins > 1 else "Variable"

    # Create histogram
    plt.figure(figsize=(10, 6))
    ax = sns.histplot(cluster_sizes, bins=nbins, kde=True)

    # Plot Title and Axes labels
    #plt.title(f'Distribution of Cluster Sizes\nNumber of Bins: {num_bins}, Bin Size: {bin_size}')
    plt.title(f"Distribution of Cluster Sizes"
          f"{f' (Clusters with >{min_proteomes} Proteomes)' if min_proteomes is not None else ''}\n"
          f"Number of Bins: {num_bins}, Bin Size: {bin_size}")

    
    plt.xlabel(f"{xlabel} {'(%)' if x_proteomes == 'percent' else ''}") 
    
    if ylog:
        plt.ylabel(f'{ylabel} (Log)')
        plt.yscale('log')
    else:
        plt.ylabel(f'{ylabel}')

    plt.tight_layout()

    # Compute and display Core/Accessory/Rare stats
    if show_stats:
        if total_proteomes is None:
            print("Warning: total_proteomes is required for core/accessory/rare classification.")
        else:
            # Convert cluster sizes to percentage presence
            percentages = (df_filtered[plot_colname] / total_proteomes) * 100

            # Compute counts for each category
            num_clusters = len(percentages)
            hard_core = (percentages == 100).sum()
            soft_core = ((percentages >= 95) & (percentages < 100)).sum()
            accessory = ((percentages >= 5) & (percentages < 95)).sum()
            rare = (percentages < 5).sum()
            
            # If filtering was applied, report original & filtered counts
            original_clusters = len(df)  # Before filtering
            min_proteomes_text = f"(Filtered: min_proteomes > {min_proteomes})" if min_proteomes is not None else ""

            # Print stats
            print("\nDescriptive Statistics:")
            print(f"  Total proteomes: {total_proteomes}")
            print(f"  Total proteins: {total_proteins}")
            #print(f"  Total Clusters: {num_clusters}")
            print(f"  Original Total Clusters: {original_clusters:,}")
            if min_proteomes is not None:
                print(f"  Clusters After Filtering: {num_clusters:,} {min_proteomes_text}")
            print(f"  Hard Core Clusters (100%): {hard_core}")
            print(f"  Soft Core Clusters (≥95%): {soft_core}")
            print(f"  Accessory Clusters (<95% & ≥5%): {accessory}")
            print(f"  Rare Clusters (<5%): {rare}")

            # Add to plot
            stats_text = (
                f"Total Proteomes: {total_proteomes:,}\n"
                f"Total Proteins: {total_proteins:,}\n"
                f"-----------------------------------\n"
                f"Total Clusters: {num_clusters:,}\n"
                f"Hard Core (100%): {hard_core:,}\n"
                f"Soft Core (≥95%): {soft_core:,}\n"
                f"Accessory (≥5%, <95%): {accessory:,}\n"
                f"Rare (<5%): {rare:,}"
            )
            
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.5, 0.8, 
                    stats_text, 
                    transform=ax.transAxes, 
                    fontsize=12,
                    verticalalignment='center', 
                    horizontalalignment='left', 
                    bbox=props)
            
    # Save or show plot
    if output_plot:
        plt.savefig(output_plot)
        print(f"\nCluster size distribution plot saved to: {output_plot}")
    else:
        plt.show()

###############################################################################
# Example usage:
#homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
#basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"
#input_file = basedir+"/DEL/clusterizes_proteomes2"
#output_plot_file = 'out2.png'

#input_df = pd.read_csv(input_file, sep='\t', names=['proteome_count'] )
#tot_proteomes = 146237
#tot_proteins = 307483396

# call function
#plot_clustersize(df=input_df, 
#          plot_colname='proteome_count', 
#          total_proteomes=tot_proteomes,
#          total_proteins=tot_proteins,
#          min_proteomes=None, #or an int like 1
#          ylog=True, 
#          x_proteomes='number', # percent
#          output_plot='out.png',
#          show_stats=True)
