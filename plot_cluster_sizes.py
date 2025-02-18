#!/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm

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
    print(f"Cluster count saved to: {output_file}")

def plot_data(input_file, output_plot, min_proteomes):
    """Plots the distribution of cluster sizes, filtering by minimum proteomes if needed."""
    df = pd.read_csv(input_file, sep='\t')

    # Filter based on min_proteomes threshold
    df = df[df['count_of_proteomes'] >= min_proteomes]
    
    plt.figure(figsize=(10, 6))
    plt.hist(df['count_of_proteomes'], bins=30, edgecolor='black')
    plt.title('Distribution of Cluster Sizes')
    plt.xlabel('Cluster Size (Number of Proteomes)')
    plt.ylabel('Frequency')
    plt.grid(True)
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

    args = parser.parse_args()

    # Set default output filenames
    if args.output_file is None:
        args.output_file = f"Count_{args.input_file}"
    if args.output_plot is None:
        args.output_plot = f"Plot_{args.input_file}.png"

    # Step 1: Process the file if requested
    if args.process:
        print(f"Processing file: {args.input_file}")
        process_file(args.input_file, args.output_file, args.chunksize)

    # Step 2: Generate plot if requested
    if args.plot:
        plot_input = args.output_file if args.process else args.input_file
        print(f"Generating plot from: {plot_input}")
        plot_data(plot_input, args.output_plot, args.min_proteomes)
