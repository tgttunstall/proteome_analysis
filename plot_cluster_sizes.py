#!/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm

input_file="/home/tunstall/Documents/arise/spneumo_dataset/outL.tsv"
chunksize=100000
#output_file="/home/tunstall/Documents/arise/spneumo_dataset/Count_outL.tsv"
output_file="/home/tunstall/Documents/arise/spneumo_dataset/Count_Labelled_Species_protein_cluster.tsv"
min_proteomes=0

#########
df = pd.read_csv(input_file, sep = "\t")
a = df.groupby(['cluster_id',  'protein_id']).size()
b = df.groupby('cluster_id')['protein_id'].nunique()
c = df.groupby('protein_id')['cluster_id'].nunique()

#########



def process_file(input_file, output_file, chunksize):
    """Processes the input file to count proteomes per cluster and saves the result."""
    cluster_counts = pd.Series(dtype=int)
    
    # Count total lines for progress tracking
    total_lines = sum(1 for _ in open(input_file, 'r'))
    progress_bar = tqdm(total=total_lines, desc="Processing", unit="lines")
    
    # Read in chunks and count occurrences
    for chunk in pd.read_csv(input_file, sep='\t', usecols=['cluster_id'], chunksize=chunksize):
        cluster_counts = cluster_counts.add(chunk['cluster_id'].value_counts(), fill_value=0)
        progress_bar.update(len(chunk))
    
    progress_bar.close()
    
    # Save to file
    cluster_counts = cluster_counts.astype(int).reset_index()
    cluster_counts.columns = ['cluster_id', 'count_of_proteomes']
    cluster_counts.to_csv(output_file, sep='\t', index=False)
    print(f"\nCluster count saved to: {output_file}")

def plot_data(input_file, output_plot, min_proteomes=0):
    """Plots the distribution of cluster sizes, filtering by minimum proteomes if needed."""
    df = pd.read_csv(input_file, sep='\t')
    #df = pd.read_csv(output_file, sep='\t')

    # Filter based on min_proteomes threshold
    df = df[df['count_of_proteomes'] >= min_proteomes]
    
    plt.figure(figsize=(10, 6))
    plt.hist(df['count_of_proteomes'], bins=30, edgecolor='black')
    plt.title('Distribution of Cluster Sizes')
    plt.xlabel('Cluster Size (Number of Proteomes)')
    plt.ylabel('Frequency')
    plt.grid(True)
    
    # Apply log scale to the axis
    #plt.yscale("log")
    #plt.xscale("log")
    
    plt.savefig(output_plot)
    print(f"Plot saved to: {output_plot}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process cluster files and generate plots.")
    parser.add_argument("input_file", help="Path to the input file (raw dataset or precomputed count file).")
    parser.add_argument("--process", action="store_true", help="Compute cluster proteome counts.")
    parser.add_argument("--plot", action="store_true", help="Generate cluster size distribution plot.")
    parser.add_argument("--output_file", default=None, help="Custom output file for counts (default: 'Count_<input_file>').")
    parser.add_argument("--output_plot", default=None, help="Custom output file for plot (default: 'Plot_<input_file>.png').")
    parser.add_argument("--chunksize", type=int, default=100000, help="Number of rows to process per chunk (default: 100000).")
    parser.add_argument("--min_proteomes", type=int, default=1, help="Minimum proteomes per cluster to include in plot (default: 1).")

    # Step 2: Generate plot if requested
    if args.plot:
        plot_input = args.output_file if args.process else args.input_file
        print(f"Generating plot from: {plot_input}")
        plot_data(plot_input, args.output_plot, args.min_proteomes)

###############################################################################:wq

###############################################################################
#input_file = "/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Labelled_Species_protein_cluster.tsv"
#output_file = "/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Count_Labelled_Species_protein_cluster.tsv"
#chunksize= 100000
#min_proteomes=0

#f = pd.read_csv("/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Labelled_Species_protein_cluster.tsv", sep = "\t")

#f.groupby(['cluster_id']).size()
#f.groupby(['protein_id']).size()
#f.groupby(['proteomes']).size()
#f.groupby(['cluster_id', 'proteomes'])
#f.value_counts(subset = ['cluster_id'])

#f2 = f['cluster_id'].value_counts()

