#!/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

def process_file(input_file, output_file):
    # Initialize an empty series to store cluster counts
    cluster_counts = pd.Series(dtype=int)

    # Determine the total number of lines to set up the progress bar accurately
    total_lines = sum(1 for line in open(input_file, 'r'))
    progress_bar = tqdm(total=total_lines, desc="Processing", unit="lines")

    # Read the file in chunks of 100,000 rows
    for chunk in pd.read_csv(input_file, sep='\t', chunksize=100000, usecols=['cluster_id']):
        cluster_counts = cluster_counts.add(chunk['cluster_id'].value_counts(), fill_value=0)
        progress_bar.update(len(chunk))

    progress_bar.close()

    # Plotting the cluster sizes as a histogram
    plt.figure(figsize=(10, 6))
    cluster_counts.hist(bins=30)
    plt.title('Distribution of Cluster Sizes')
    plt.xlabel('Cluster Size')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(output_file)
    print(f"Graph has been saved to {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.tsv output_graph.png")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        print(f"Starting to process the file: {input_file}")
        process_file(input_file, output_file)
        print(f"Processing complete. Graph has been saved as: {output_file}")


input_file = "/home/tunstall/Documents/arise/spneumo_dataset/outL.tsv"
