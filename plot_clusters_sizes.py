#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from tqdm import tqdm
import argparse

def process_file(input_file, output_file, chunksize):
    """
    Reads a large input file in chunks, counts unique proteomes per cluster_id, 
    and writes the counts to the output file.
    """
    # Initialize a dictionary to store cluster_id → unique proteome count
    cluster_proteome_count = {}

    # Determine the total number of lines for progress bar
    with open(input_file, 'r') as f:
        total_lines = sum(1 for _ in f)

    progress_bar = tqdm(total=total_lines, desc="Processing", unit="lines")

    # Read file in chunks, ensuring headers are correctly used
    col_names = ['cluster_id', 'protein_id', 'proteomes', 'is_rep']
    
    for chunk in pd.read_csv(input_file, sep='\t', usecols=['cluster_id', 'proteomes'], header=0, chunksize=chunksize):
        # Count unique proteomes per cluster_id
        cluster_counts = chunk.groupby('cluster_id')['proteomes'].nunique()

        # Merge with existing counts
        for cluster, count in cluster_counts.items():
            cluster_proteome_count[cluster] = cluster_proteome_count.get(cluster, 0) + count

        progress_bar.update(len(chunk))

    progress_bar.close()

    # Convert results to DataFrame and save
    count_df = pd.DataFrame(list(cluster_proteome_count.items()), columns=['cluster_id', 'count_of_proteomes'])
    count_df.sort_values(by='cluster_id', inplace=True)
    count_df.to_csv(output_file, sep='\t', index=False)

    print(f"Count data saved to: {output_file}")

def plot_data(csv_file, x_col, y_col, output_plot):
    """
    Plots data from the given CSV file using specified x and y columns.
    """
    df = pd.read_csv(csv_file, sep='\t')
    
    if x_col not in df.columns or y_col not in df.columns:
        print(f"Error: Specified columns '{x_col}' and '{y_col}' not found in {csv_file}.")
        print(f"Available columns: {df.columns.tolist()}")
        return

    plt.figure(figsize=(10, 6))
    plt.bar(df[x_col], df[y_col])
    plt.xlabel(x_col.replace("_", " ").title())
    plt.ylabel(y_col.replace("_", " ").title())
    plt.title(f'{y_col.replace("_", " ").title()} per {x_col.replace("_", " ").title()}')
    plt.grid(True)
    plt.savefig(output_plot)
    print(f"Plot saved to: {output_plot}")

if __name__ == '__main__':
    # Argument parsing
    parser = argparse.ArgumentParser(description="Count unique proteomes per cluster and optionally plot results.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("-o", "--output_file", help="Path to output file (default: Count_<input_file>)", default=None)
    parser.add_argument("--plot", action="store_true", help="Generate a plot (default: False)")
    parser.add_argument("--output_plot", help="Path to save the plot (default: Plot_<input_file>.png)", default=None)
    parser.add_argument("--chunksize", type=int, default=100000, help="Number of rows per chunk (default: 100000)")
    parser.add_argument("--x_col", default="cluster_id", help="Column name for x-axis (default: cluster_id)")
    parser.add_argument("--y_col", default="count_of_proteomes", help="Column name for y-axis (default: count_of_proteomes)")

    args = parser.parse_args()

    # Set default output filename if not provided
    if args.output_file is None:
        args.output_file = f"Count_{os.path.basename(args.input_file)}"

    print(f"Processing file: {args.input_file} with chunksize={args.chunksize}")
    
    # Step 1: Process the file to count unique proteomes per cluster
    process_file(args.input_file, args.output_file, args.chunksize)

    # Step 2: If plotting is enabled, call plot function with user-specified columns
    if args.plot:
        output_plot = args.output_plot or f"Plot_{os.path.basename(args.input_file)}.png"
        plot_data(args.output_file, args.x_col, args.y_col, output_plot)

    print("Processing complete.")
