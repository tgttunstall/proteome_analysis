#!/usr/bin/env python3
import time
import os
import argparse
from glob import glob
import pprint as pp

# Set up argument parser
parser = argparse.ArgumentParser(description="Proteome comparison script with configurable input and output paths.")
parser.add_argument("--fasta_dir", type=str, required=True, help="Directory containing all .fa files")
parser.add_argument("--tsv_file", type=str, required=True, help="Path to the TSV file containing protein clusters")
parser.add_argument("--out_file", type=str, required=True, help="Path to the output file")

# Parse arguments
args = parser.parse_args()
proteome_indir = args.fasta_dir
protein_cluster_infile = args.tsv_file
output_file = args.out_file

#proteome_indir = "/home/tunstall/Documents/arise/data_arise_proteome/pig/"
#protein_cluster_infile = "/home/tunstall/Documents/arise/data_arise_proteome/results_pig/Specie_protein_cluster.tsv"
#output_file = "/home/tunstall/Documents/arise/data_arise_proteome/results_pig/output.tsv"


# Start timing
start_time = time.time()

#===============
# Stage 1: Create proteome dictionary
#===============

# Initialise the protein-to-proteome dictionary
proteomeD = {}

# Find all .fa files
fa_files = glob(os.path.join(proteome_indir, "*.fa"))
print(f"Number of .fa files found: {len(fa_files)}")
if len(fa_files) == 0:
    print("No .fa files found. Please check the directory path.")

# Populate proteomeD directly from .fa files
for file in fa_files:
    proteome_id = os.path.splitext(os.path.basename(file))[0]

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line.strip().lstrip('>')
                if protein_id not in proteomeD:
                    proteomeD[protein_id] = proteome_id
                else:
                    proteomeD[protein_id] += f"_{proteome_id}"

print("Sample dictionary entries:")
pp.pprint(list(proteomeD.items())[:5])

# End timing for stage 1
end_time = time.time()
elapsed_time_s1 = end_time - start_time
print(f"Dictionary creation completed in {elapsed_time_s1 // 3600:.0f} hours, "
      f"{(elapsed_time_s1 % 3600) // 60:.0f} minutes, and {elapsed_time_s1 % 60:.2f} seconds.")

#===============
# Stage 2: Assign proteome ID to protein names
#===============

# Helper function to look up the proteome ID
def find_proteome_id(protein_name):
    return proteomeD.get(protein_name, "XXX")

# Read in the TSV file with protein pairs
with open(protein_cluster_infile, 'r') as tsvfile:
    unique_entries = set()
    for line in tsvfile:
        cols = line.strip().split('\t')
        if len(cols) >= 2:
            protein1_id, protein2_id = cols[0], cols[1]
            unique_entries.add((protein1_id, protein2_id))

print(f"\nLength of unique entries: {len(unique_entries)}")

# Process the unique entries
cluster_data = []
cluster_id = 0
seen_clusters = {}

for protein1_id, protein2_id in unique_entries:
    if protein1_id not in seen_clusters:
        cluster_id += 1
        seen_clusters[protein1_id] = cluster_id
        representative = '*'
    else:
        cluster_id = seen_clusters[protein1_id]
        representative = ''

    proteome_id_2 = find_proteome_id(protein2_id)
    protein2_id_proteome = f"{proteome_id_2}_{protein2_id}"
    cluster_data.append((cluster_id, representative, protein2_id_proteome))

# End timing for stage 2
end_time = time.time()
elapsed_time_s2 = end_time - start_time
print(f"Mapping completed in {elapsed_time_s2 // 3600:.0f} hours, "
      f"{(elapsed_time_s2 % 3600) // 60:.0f} minutes, and {elapsed_time_s2 % 60:.2f} seconds.")

#===============
# Write Output File
#===============

print(f"Writing file: {output_file}\nNumber of clusters: {len(cluster_data)}")
with open(output_file, 'w') as outfile:
    outfile.write("cluster_id\trepresentative\tprotein2_id_proteome\n")
    for row in cluster_data:
        outfile.write('\t'.join(map(str, row)) + '\n')

print(f"Updated file saved as {output_file}")

