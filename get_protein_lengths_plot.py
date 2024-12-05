#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 08:55:36 2024

@author: tanu
"""

#from Bio import SeqIO

# # Path to your multi-line FASTA file
# fasta_file = "/home/pub/Work/data_arise_proteome/ggcaller/gene_calls.faa"

# # Read the file and calculate sequence lengths
# for record in SeqIO.parse(fasta_file, "fasta"):
#     seq_id = record.id  # Get the sequence ID
#     seq_length = len(record.seq)  # Get the length of the sequence
#     print(f"Sequence ID: {seq_id}, Length: {seq_length}")
    
########
import time
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#fasta_file = "/home/pub/Work/data_arise_proteome/ggcaller/gene_calls.faa"
fasta_file = "/home/pub/Work/data_arise_proteome/ggcaller/cleaned_gene_calls.faa"

output_file = "/home/pub/Work/data_arise_proteome/ggcaller/BP_sequence_lengths2.txt"

# Start timing
start_time = time.time()

# # Initialize variables
# zero_length_sequences = []
# total_length = 0
# num_sequences = 0
# min_length = float('inf')  # Start with a very high value
# max_length = 0  # Start with a very low value

# # Count total sequences for tqdm progress bar
# total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

# # Open output file
# with open(output_file, "w") as out_f:
#     out_f.write("Header\tLength\n")  # Write header for the output file
    
#     # Process the FASTA file with a progress bar
#     for record in tqdm(SeqIO.parse(fasta_file, "fasta"), total=total_sequences, desc="Processing sequences"):
#         seq_id = record.id
#         seq_length = len(record.seq)
        
#         # Alert for zero-length sequences
#         if seq_length == 0:
#             zero_length_sequences.append(seq_id)
        
#         # Update min and max
#         if seq_length < min_length:
#             min_length = seq_length
#         if seq_length > max_length:
#             max_length = seq_length
        
#         # Write to output file
#         out_f.write(f"{seq_id}\t{seq_length}\n")
        
#         # Update statistics
#         total_length += seq_length
#         num_sequences += 1

# # Calculate summary statistics
# end_time = time.time()
# average_length = total_length / num_sequences if num_sequences > 0 else 0

# # Print results
# print(f"Processed {num_sequences} sequences.")
# print(f"Minimum sequence length: {min_length}")
# print(f"Maximum sequence length: {max_length}")
# print(f"Average sequence length: {average_length:.2f}")
# print(f"Total time: {end_time - start_time:.2f} seconds.")

# if zero_length_sequences:
#     print(f"Found {len(zero_length_sequences)} zero-length sequences. Headers:")
#     for header in zero_length_sequences:
#         print(f" - {header}")
# else:
#     print("No zero-length sequences found.")


def calculate_fasta_length(fasta_file, output_file, total_sequences=None):
    """
    Process a multi-line FASTA file to calculate sequence lengths.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_file (str): Path to save the output lengths.
        total_sequences (int, optional): Total number of sequences for progress tracking. Defaults to None.

    Returns:
        dict: Summary statistics (total, average, min, max, zero_length_headers).
    """
    start_time = time.time()

    # Initialize variables
    zero_length_sequences = []
    total_length = 0
    num_sequences = 0
    min_length = float('inf')  # Start with a very high value
    max_length = 0  # Start with a very low value

    # Calculate total sequences if not supplied
    if total_sequences is None:
        print("Calculating total sequences...")
        total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        print(f"Total sequences: {total_sequences}")

    # Open output file
    with open(output_file, "w") as out_f:
        out_f.write("Header\tLength\n")  # Write header for the output file

        # Process the FASTA file with a progress bar
        for record in tqdm(SeqIO.parse(fasta_file, "fasta"), total=total_sequences, desc="Processing sequences"):
            seq_id = record.id
            seq_length = len(record.seq)

            # Alert for zero-length sequences
            if seq_length == 0:
                zero_length_sequences.append(seq_id)

            # Update min and max
            if seq_length < min_length:
                min_length = seq_length
            if seq_length > max_length:
                max_length = seq_length

            # Write to output file
            out_f.write(f"{seq_id}\t{seq_length}\n")

            # Update statistics
            total_length += seq_length
            num_sequences += 1

    # Calculate summary statistics
    end_time = time.time()
    average_length = total_length / num_sequences if num_sequences > 0 else 0

    # Print results
    print(f"\nProcessed {num_sequences} sequences.")
    print(f"Minimum sequence length: {min_length}")
    print(f"Maximum sequence length: {max_length}")
    print(f"Average sequence length: {average_length:.2f}")
    print(f"Total time: {end_time - start_time:.2f} seconds.")

    if zero_length_sequences:
        print(f"Found {len(zero_length_sequences)} zero-length sequences.")
    else:
        print("No zero-length sequences found.")

    # Return summary statistics
    return {
        "total_sequences": num_sequences,
        "average_length": average_length,
        "min_length": min_length,
        "max_length": max_length,
        "zero_length_headers": zero_length_sequences,
        "processing_time": end_time - start_time,
    }




#summary_stats = calculate_fasta_length(fasta_file = fasta_file, 
#                                       output_file = output_file,
#                                       total_sequences = 9496341)


def process_fasta_with_plot(fasta_file, output_file, total_sequences=None, plot=False):
    """
    Process a multi-line FASTA file and optionally plot a histogram of sequence lengths.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_file (str): Path to save the output lengths.
        total_sequences (int, optional): Total number of sequences for progress tracking. Defaults to None.
        plot (bool): Whether to plot a histogram of sequence lengths. Defaults to False.

    Returns:
        dict: Summary statistics (total, average, min, max, SD, zero_length_headers).
    """
    start_time = time.time()

    # Initialize variables
    zero_length_sequences = []
    total_length = 0
    num_sequences = 0
    min_length = float('inf')  # Start with a very high value
    max_length = 0  # Start with a very low value
    lengths = []  # Collect all sequence lengths

    # Calculate total sequences if not supplied
    if total_sequences is None:
        print("Calculating total sequences...")
        total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        print(f"Total sequences: {total_sequences}")

    # Open output file
    with open(output_file, "w") as out_f:
        out_f.write("Header\tLength\n")  # Write header for the output file

        # Process the FASTA file with a progress bar
        for record in tqdm(SeqIO.parse(fasta_file, "fasta"), total=total_sequences, desc="Processing sequences"):
            seq_id = record.id
            seq_length = len(record.seq)
            lengths.append(seq_length)

            # Alert for zero-length sequences
            if seq_length == 0:
                zero_length_sequences.append(seq_id)

            # Update min and max
            if seq_length < min_length:
                min_length = seq_length
            if seq_length > max_length:
                max_length = seq_length

            # Write to output file
            out_f.write(f"{seq_id}\t{seq_length}\n")

            # Update statistics
            total_length += seq_length
            num_sequences += 1

    # Calculate summary statistics
    end_time = time.time()
    average_length = total_length / num_sequences if num_sequences > 0 else 0
    std_dev = np.std(lengths) if num_sequences > 0 else 0

    # Print results
    print(f"\nProcessed {num_sequences} sequences.")
    print(f"Minimum sequence length: {min_length}")
    print(f"Maximum sequence length: {max_length}")
    print(f"Average sequence length: {average_length:.2f}")
    print(f"Standard deviation: {std_dev:.2f}")
    print(f"Total time: {end_time - start_time:.2f} seconds.")

    if zero_length_sequences:
        print(f"Found {len(zero_length_sequences)} zero-length sequences.")
    else:
        print("No zero-length sequences found.")

    # Plot histogram if requested
    if plot:
        # Plot histogram of sequence lengths
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, edgecolor="black")
        plt.title("Histogram of Sequence Lengths")
        plt.xlabel("Sequence Length")
        plt.ylabel("Frequency")

        # Annotate summary statistics
        stats_text = (
            f"Min: {min_length}\n"
            f"Max: {max_length}\n"
            f"Avg: {average_length:.2f}\n"
            f"SD: {std_dev:.2f}\n"
            f"Total: {total_sequences}"
        )
        plt.text(0.95, 0.95, stats_text, fontsize=10, transform=plt.gca().transAxes,
                 verticalalignment='top', horizontalalignment='right',
                 bbox=dict(facecolor='white', alpha=0.5))

        plt.grid(axis="y")
        plt.show()

    # Return summary statistics
    return {
        "total_sequences": num_sequences,
        "average_length": average_length,
        "min_length": min_length,
        "max_length": max_length,
        "standard_deviation": std_dev,
        "zero_length_headers": zero_length_sequences,
        "processing_time": end_time - start_time,
    }

#total_sequences = 9496341
summary_stats = process_fasta_with_plot(fasta_file, output_file, total_sequences=None, plot=True)
