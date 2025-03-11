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
              output_plot=None, 
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
###############################################################################
def plot_category_counts(df,
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
                         plot_title="Distribution of Total and Exclusive Counts",
                         plot_label_size=16,
                         legend_font_size=12,
                         output_plot=None):
    """
    Generates a barplot of cluster sizes from a precomputed counts file.

    Args:
        df: DataFrame containing the data.
        
    """
    colors = cat_colours
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot total counts
    bars_total = ax.bar(df[category_col], df[total_col], color=[colors[cat] for cat in df[category_col]], label=total_col_label)

    # Plot exclusive counts with hatching
    bars_exclusive = ax.bar(df[category_col], df[exclusive_col], color=[colors[cat] for cat in df[category_col]],
                            hatch="///", edgecolor="black", alpha=1, label=exclusive_col_label)

    # Annotate values on the bars
    for bar, value in zip(bars_total, df[total_col]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1000, f"{value:,}", 
                ha="center", fontsize=bar_annot_size, color="black")

    for bar, value in zip(bars_exclusive, df[exclusive_col]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() - 10000, f"{value:,}", 
                ha="center", fontsize=bar_annot_size, color="white", fontweight="bold")

    # Set labels and title
    ax.set_xlabel(xlabel, fontsize=label_size)
    ax.set_ylabel(ylabel, fontsize=label_size)
    ax.set_title(plot_title, fontsize=plot_label_size)

    # Custom legend to indicate Total and Exclusive
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor="grey", edgecolor="black", label=total_col_label),
        plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="black", hatch="///", label=exclusive_col_label)
    ]
    ax.legend(handles=legend_elements, fontsize=legend_font_size)
    # Save or show plot
    if output_plot:
        plt.savefig(output_plot)
        print(f"\nCluster size distribution plot saved to: {output_plot}")
    else:
        plt.show()
    
# Example usage:
#aC = 26000
#bC = 50000
#cC = aC+bC
#eg_data = {'Categories': ['C1', 'C2', 'C3'],
#        'Total': [aC, bC, cC],
#        'Exclusive': [1000, 5000, 8000]}

# Create DataFrame
#bar_df = pd.DataFrame(eg_data)

# print(bar_df)
# plot_category_counts(df = bar_df,
#                      category_col="Categories", 
#                      total_col="Total", 
#                      total_col_label = "Total",
#                      exclusive_col="Exclusive",
#                      exclusive_col_label = "Exclusive",
#                      cat_colours = {"C1": "blue",
#                                    "C2": "darkorange",
#                                    "C3": "purple"},
#                      bar_annot_size = 12,
#                      xlabel = "",
#                      ylabel = "Count",
#                      label_size = 14,
#                      plot_title="Example",
#                      plot_label_size=16,
#                      legend_font_size=12,
#                      output_plot=False)

###############################################################################
#!/usr/bin env python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#########
def plot_histogram(data,
                   x,
                   #kde=True,
                   #nbins=None,
                   bin_method='rice',
                   color='grey',
                   label='T', 
                   alpha=0.7,
                   ylog=True,
                   xlabel='x',
                   ylabel='Count',
                   extra_title_label='',
                   output_plot=None):
    """
    Generates a histogram using Seaborn, with automatic bin calculation and NaN handling.

    Parameters
    ----------
    data : DataFrame
        Input dataset.
    x : str
        Column name for the x-axis.
    bin_method : str, optional
        Method for calculating bins (default is 'rice').
    color : str, optional
        Color of the histogram bars (default is 'grey').
    label : str, optional
        Label for legend (default is 'T').
    alpha : float, optional
        Transparency level for bars (default is 0.7).
    ylog : bool, optional
        Whether to use logarithmic scale on the y-axis (default is True).
    xlabel : str, optional
        X-axis label (default is 'x').
    ylabel : str, optional
        Y-axis label (default is 'Count').
    extra_title_label : str, optional
        Extra text to add to the plot title (default is empty).
    output_plot : str, optional
        If provided, saves the plot to this filename.

    Returns
    -------
    None
    """

    # Handle bin calculation
    # if nbins is not None:
    #     print(f"\nUsing specified number of bins: {nbins}. Ignoring bin_method '{bin_method}'.")

    #     # Calculate bin edges and bin width dynamically
    #     data_min, data_max = data[x].min(), data[x].max()
    #     bins = np.linspace(data_min, data_max, nbins + 1)  # Evenly spaced bins
    #     bin_size = round((data_max - data_min) / nbins, 2) if nbins > 1 else "Variable"

    # else:
    #     bins = np.histogram_bin_edges(data[x], bins=bin_method)
    #     print(f"\nCalculating bin edges using method: '{bin_method}'.")

    #     num_bins = len(bins) - 1
    #     bin_size = round(bins[1] - bins[0], 2) if num_bins > 1 else "Variable"
        

    # Handle missing values (omit NaNs in column `x`)
    initial_count = len(data)
    data = data.dropna(subset=[x])
    final_count = len(data)
    removed_nans = initial_count - final_count

    if removed_nans > 0:
        print(f"⚠️ Alert: Removed {removed_nans} NaN values from '{x}' column.")

    # Calculate bin edges
    bins = np.histogram_bin_edges(data[x], bins=bin_method)
    print(f"\nCalculating bin edges using method: '{bin_method}'.")
    
    num_bins = len(bins) - 1
    bin_size = round(bins[1] - bins[0], 2) if num_bins > 1 else "Variable"

    # Create histogram plot
    plt.figure(figsize=(10, 6))
    ax = sns.histplot(
        data=data,
        x=x,
        bins=num_bins,
        color=color,
        label=label,
        alpha=alpha
    )

    # Set plot titles and labels
    plt.title(f"Distribution of {xlabel}\n{extra_title_label} (n={final_count:,}) | Bins: {num_bins}, Bin Size: {bin_size}")
    plt.xlabel(xlabel)
    
    if ylog:
        plt.ylabel(f'{ylabel} (Log)')
        plt.yscale('log')
    else:
        plt.ylabel(ylabel)
    
    plt.tight_layout()

    # Save or show the plot
    if output_plot:
        plt.savefig(output_plot)
        print(f"\nPlot of Distribution of {xlabel} saved to: {output_plot}")
    else:
        plt.show()
###############################################################################
def plot_protein_busco(df, 
                       spacing_factor=1, 
                       bar_width= 0.15,
                       gap= 0.18,
                       title_fontsize=14, 
                       label_fontsize=14, 
                       tick_fontsize=12, 
                       legend_fontsize=14, 
                       fig_height = 10,
                       plot_title="Protein count and BUSCO: Duplicate samples",
                       output_plot=None):
    """
     Create a bar plot comparing protein counts and BUSCO scores for duplicate biosamples.
    
     Parameters:
     -----------
     df : pandas.DataFrame
         Input DataFrame containing biosample data with columns:
         'biosample', 'protein_count1', 'protein_count2', 'complete_combined_score1', 
         'complete_combined_score2', 'is_excluded1', 'is_excluded2', 'is_effective1', 'is_effective2'.
    
     spacing_factor : float, optional (default=1)
         Factor to adjust spacing between biosample groups on the x-axis.
    
     bar_width : float, optional (default=0.15)
         Width of each bar in the plot.
    
     gap : float, optional (default=0.18)
         Gap between bars within each biosample group.
    
     title_fontsize : int, optional (default=14)
         Font size for the plot title.
    
     label_fontsize : int, optional (default=14)
         Font size for axis labels.
    
     tick_fontsize : int, optional (default=12)
         Font size for tick labels.
    
     legend_fontsize : int, optional (default=14)
         Font size for legend (if used).
    
     fig_height : int, optional (default=10)
         Height of the figure in inches.
    
     plot_title : str, optional (default="Protein Count and BUSCO score for duplicate Biosamples")
         Title of the plot.
    
     Returns:
     --------
     fig : matplotlib.figure.Figure
         The created figure object.
    
     ax1 : matplotlib.axes.Axes
         The primary y-axis for protein counts.
    
     ax2 : matplotlib.axes.Axes
         The secondary y-axis for BUSCO scores.
    
     Notes:
     ------
     - The function creates a grouped bar plot with protein counts on the left y-axis and BUSCO scores on the right y-axis.
     - Bars are colored based on 'is_excluded' and 'is_effective' values:
       - Green: 'is_excluded' is 'f' and 'is_effective' is 'False'
       - Red: 'is_excluded' is 't' and 'is_effective' is 'True'
       - Gold: Otherwise
     - Bars for the second set of data (columns ending with '2') are hatched for differentiation.
     - Text labels are added above each bar showing the exact values.
     - X-axis labels (biosamples) are rotated 45 degrees for better readability.
    
     Example:
     --------
     df : pandas.DataFrame
    Input DataFrame containing biosample data with the following columns:
    - 'biosample': Unique identifier for each biosample.
    - 'upid1', 'upid2': IDs for duplicate samples.
    - 'protein_count1', 'protein_count2': Protein counts for each duplicate biosample.
    - 'complete_combined_score1', 'complete_combined_score2': BUSCO scores for each duplicate biosample.
    - 'is_excluded1', 'is_excluded2': Flags indicating exclusion status ('t' or 'f').
    - 'is_effective1', 'is_effective2': Flags indicating effectiveness ('True' or 'False').

    Example Data Format:
    --------------------
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
        "is_effective1": ["True", "False"],
        "is_effective2": ["False", "False"],
    }

     >>> fig, ax1, ax2 = plot_protein_busco(df, spacing_factor=1.2, bar_width=0.2, gap=0.2, title_fontsize=16)
     >>> plt.show()
     """

    # Calculate figure size based on number of samples
    fig_width = max(18, len(df) * 2)  # Minimum width of 18, scales with data
    fig_height = fig_height

    # Define positions and bar width
    x = np.arange(len(df)) * spacing_factor  # Adjust spacing between biosamples
    bar_width = bar_width
    gap = gap  # Slight gap between '1' and '2' groups

# FIXME: do this outside
# FIXME: as the range param interefers with the plot if indices are not from 0..9:
    # Function to determine bar color based on conditions
    def get_bar_color(is_excluded, is_effective):
        if is_excluded == "f" and is_effective == 'False':
            return "green"
        elif is_excluded == "t" and is_effective == 'True':
            return "red"
        else:
            return "gold"

    # Create the figure and axis objects
    fig, ax1 = plt.subplots(figsize=(fig_width, fig_height))
    ax2 = ax1.twinx()

    # Plot protein_count1 and complete_combined_score1 for each biosample (left group)
    ax1.bar(
        x - bar_width - gap / 2,
        df["protein_count1"],
        bar_width,
        edgecolor='blue',
        label="Protein Count 1",
        color=[get_bar_color(df["is_excluded1"][i], df["is_effective1"][i]) for i in range(len(df))],
    )
    ax2.bar(
        x - gap / 2,
        df["complete_combined_score1"],
        bar_width,
        edgecolor='black',
        label="BUSCO Score 1",
        color=[get_bar_color(df["is_excluded1"][i], df["is_effective1"][i]) for i in range(len(df))],
    )

    # Plot protein_count2 and complete_combined_score2 for each biosample (right group)
    ax1.bar(
        x + gap / 2,
        df["protein_count2"],
        bar_width,
        label="Protein Count 2",
        edgecolor='blue',
        color=[get_bar_color(df["is_excluded2"][i], df["is_effective2"][i]) for i in range(len(df))],
        hatch="//",  # Add hatches to differentiate columns ending in '2'
    )
    ax2.bar(
        x + bar_width + gap / 2,
        df["complete_combined_score2"],
        bar_width,
        label="BUSCO Score 2",
        edgecolor='black',
        color=[get_bar_color(df["is_excluded2"][i], df["is_effective2"][i]) for i in range(len(df))],
        hatch="//",  # Add hatches to differentiate columns ending in '2'
    )

    # Set labels for the axes with custom font sizes
    ax1.set_ylabel("Protein Count", color="blue", fontsize=label_fontsize)
    ax2.set_ylabel("BUSCO Score", color="black", fontsize=label_fontsize)

    # Customize tick parameters with font size
    ax1.tick_params(axis="y", labelcolor="blue", labelsize=tick_fontsize)
    ax2.tick_params(axis="y", labelcolor="black", labelsize=tick_fontsize)
    ax1.tick_params(axis="x", labelsize=tick_fontsize)

    # Add text labels above bars for protein counts and BUSCO scores
    for i in range(len(df)):
        # Protein Count 1 (vertical display)
        ax1.text(
            x[i] - bar_width - gap / 2,
            df["protein_count1"][i] + 10,
            str(df["protein_count1"][i]),
            ha="center",
            va="bottom",
            color="blue",
            rotation=90,  # Vertical text
            fontsize=tick_fontsize,
        )
        
        # BUSCO Score 1
        ax2.text(
            x[i] - gap / 2,
            df["complete_combined_score1"][i] + 0.5,
            str(df["complete_combined_score1"][i]),
            ha="center",
            va="bottom",
            color="black",
            fontsize=tick_fontsize,
        )
        
        # Protein Count 2 (vertical display)
        ax1.text(
            x[i] + gap / 2,
            df["protein_count2"][i] + 10,
            str(df["protein_count2"][i]),
            ha="center",
            va="bottom",
            color="blue",
            rotation=90,  # Vertical text
            fontsize=tick_fontsize,
        )
        
        # BUSCO Score 2
        ax2.text(
            x[i] + bar_width + gap / 2,
            df["complete_combined_score2"][i] + 0.5,
            str(df["complete_combined_score2"][i]),
            ha="center",
            va="bottom",
            color="black",
            fontsize=tick_fontsize,
        )

    # Set x-axis labels with rotation and font size
    plt.xticks(x, df["biosample"], rotation=45, ha="right", fontsize=tick_fontsize)

    # Set plot title with custom font size
    plt.title(plot_title, fontsize=title_fontsize)

    # Adjust layout to fit rotated labels
    plt.tight_layout()
    
    # Save or show the plot
    if output_plot:
        plt.savefig(output_plot, bbox_inches="tight")
        print(f"\nPlot saved to: {output_plot}")
    else:
        plt.show()

    return fig, ax1, ax2

# Usage example:
# data = {
#     "biosample": ["SAMEA1024557", "SAMEA1024588"],
#     "upid1": ["UP000047540", "UP000235499"],
#     "upid2": ["UP000235454", "UP000046519"],
#     "protein_count1": [2010, 2028],
#     "protein_count2": [1982, 2012],
#     "complete_combined_score1": [99.8, 99.8],
#     "complete_combined_score2": [99.8, 99.3],
#     "is_excluded1": ["t", "f"],
#     "is_excluded2": ["f", "f"],
#     "is_effective1": ["True", "False"],
#     "is_effective2": ["False", "False"],
# }

# # Create DataFrame
# df = pd.DataFrame(data)
# # fig, ax1, ax2 = plot_protein_busco(df)
# # plt.show()
# fig, ax1, ax2 = plot_protein_busco(df, output_plot="protein_busco_plot.png")

