#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:16:53 2024

@author: tanu
"""
#!/usr/bin/env python3
import os
import argparse
from multiprocessing import Pool, current_process
import subprocess
import glob
import time
import sys
from typing import List, Optional
from tqdm import tqdm

#from ffdb import split_file  # Importing split_file from ffdb
from ffdb import *

################
BUFFERSIZE = 1048576  # 1Mb
DESCRIPTION = """
Script to parallelise clustal omega alignments
"""

#################
# Clustal Omega configuration
CLUSTALO_EXE = "clustalo"  # Path or name of Clustal Omega executable

# Helper Functions

def secs2time(secs):
    """Converts a time duration in seconds to a human-readable format (hours, minutes, seconds)."""
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)

def elapsed_time(start_time):
    """Computes the elapsed time from a given start time in seconds and returns a formatted string."""
    process_time = time.time() - start_time
    return secs2time(process_time)

def exit_with_error(message: str, code: int = 1):
    """Prints an error message to stderr and exits the program with the specified exit code."""
    eprint(f" => {message}")
    sys.exit(code)

def eprint(*myargs, **kwargs):
    """Prints the provided arguments to stderr."""
    print(*myargs, file=sys.stderr, **kwargs)


def output_time(start_time):
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    total_elapsed_time = elapsed_time(start_time)
    eprint(" '-- Processing complete --'")
    eprint(f" |-- END TIME: {end_time}")
    eprint(f" |-- TOTAL ELAPSED TIME: {total_elapsed_time}")    

def check_args():
    """
    Parse arguments and check for error conditions.
    """

    def positive_integer(value):
        try:
            value = int(value)
            if value <= 0:
                raise argparse.ArgumentTypeError(
                    "{} is not a positive integer".format(value)
                )
        except ValueError:
            raise argparse.ArgumentTypeError("{} is not an integer".format(value))
        return value

    def is_valid_file(path):
        """Check if the given path is a valid, readable file."""
        if not path:
            raise argparse.ArgumentTypeError(f"File path cannot be empty or None.")

        if not os.path.exists(path):
            raise argparse.ArgumentTypeError(f"The file '{path}' does not exist.")

        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"The path '{path}' is not a valid file.")

        if not os.access(path, os.R_OK):
            raise argparse.ArgumentTypeError(f"The file '{path}' is not readable.")

        return path

    class CustomArgumentParser(argparse.ArgumentParser):
        def print_help(self, *args, **kwargs):
            """
            Print custom text before the default help message.
            """
            print("Run Clustal Omega in parallel across multiple input files.")
            super().print_help(*args, **kwargs)

    parser = CustomArgumentParser(description="Parallelising Clustal Omega alignment")

    #==================================
    # General script-level arguments
    #==================================
    
    # General input options
    input_group = parser.add_argument_group("Input Options (choose one)")
    input_group.add_argument("-d", "--input_fasta_dir", type=str, required=False,
        help="Directory containing all proteome .fa files.")
    input_group.add_argument("-f", "--input_fasta_list", type=is_valid_file, required=False,
        help="Path to a text file containing a list of FASTA files.")
         
    # General options
    parser.add_argument("-t", "--threads", type=positive_integer, default=1,
        help="Number of threads for parallel processing.")
    parser.add_argument("-e", "--extension", type=str, default=".fa",
        help="File extension for FASTA files (default: .fa).")

    parser.add_argument("--sequence_stats", choices=["none", "basic", "detailed"], default="none",
                    help="Level of sequence statistics to report (none, basic, or detailed)")
    # parser.add_argument("-d", "--input_fasta_dir", type=str, required=False,
    #     help="Directory containing all proteome .fa files.")

    # parser.add_argument("-f", "--input_fasta_list", type=is_valid_file, required=False,
    #     help="Path to a text file containing a list of FASTA files.")  

    #==================================
    # Clustal Omega-specific arguments
    #==================================

    clustalo_group = parser.add_argument_group("Clustal Omega specific arguments")
    
    clustalo_group.add_argument("-o", "--out_dir", type=str, required=True,
        help="Output directory for aligned files.")

    clustalo_group.add_argument("-at", "--align_threads", type=positive_integer, default=1,
        help="Number of threads for Clustal Omega.")

    clustalo_group.add_argument("--force", action="store_true",
        help="Force overwrite of existing output files.")

    clustalo_group.add_argument("--seqtype", type=str, choices=["Protein", "DNA"], default="Protein",
        help="Specify the sequence type (e.g., 'Protein' or 'DNA').")
    
    clustalo_group.add_argument("--outfmt", type=str, choices=["fa","clu","msf","phy","selex","st","vie"], default="fa",
        help="Specify the MSA output format.")

    args = parser.parse_args()

    # Validate input arguments
    if args.input_fasta_dir and args.input_fasta_list:
        exit_with_error("ERROR: Please specify either --input_fasta_dir or --input_fasta_list, not both.", 1)

    if args.input_fasta_dir and not os.path.isdir(args.input_fasta_dir):
        exit_with_error(f"ERROR: No such directory '{args.input_fasta_dir}'", 2)

    if args.input_fasta_list and not os.path.isfile(args.input_fasta_list):
        exit_with_error(f"ERROR: No such file '{args.input_fasta_list}'", 2)

    if not os.path.isdir(args.out_dir):
        try:
            os.makedirs(args.out_dir)
            eprint(f" |-- Created output directory: {args.out_dir}")
        except PermissionError:
            exit_with_error(f"ERROR: Cannot create directory '{args.out_dir}'", 1)

    eprint(f" |-- Output directory: {args.out_dir}")
        
    #args = parser.parse_args()

    if args.threads > 1 and args.input_fasta_list:
        args.chunksize = calculate_blocksize(args.input_fasta_list, args.threads)
        eprint(
            f" |-- threads: {args.threads}; input file will be split in chunks of max {args.chunksize} bytes"
        )

    return args

def get_sequence_statistics(fasta_file):
    sequence_count = 0
    lengths = []

    with open(fasta_file, 'r') as f:
        current_length = 0
        for line in f:
            if line.startswith(">"):
                if current_length > 0:
                    sequence_count += 1
                    lengths.append(current_length)
                current_length = 0
            else:
                current_length += len(line.strip())

        if current_length > 0:
            sequence_count += 1
            lengths.append(current_length)

    return {
        "total_sequences": sequence_count,
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "avg_length": sum(lengths) / sequence_count if sequence_count > 0 else 0,
        "sequence_lengths": lengths
    }

def get_alignment(input_file, args):
    """
    Align sequences using Clustal Omega and write to the output file.
    """
    # Derive the output file path from args.out_dir
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(args.out_dir, f"{base_name}.aln")

    # Construct Clustal Omega command
    clustalo_call = [
        CLUSTALO_EXE,
        "--infile", input_file,   #f"--infile={input_file}",
        "--outfile", output_file, #"--outfile", "{}".format(output_file)
        "--outfmt", args.outfmt,
        "--threads", str(args.align_threads), #because subprocess.run expects it!
        "--seqtype", args.seqtype,
    ]

    # Add --force if specified
    if args.force:
        clustalo_call.append("--force")

    # Log the command and the number of threads being used
    #eprint(f"\nForce overwrite enabled: {args.force}")
    eprint(f"Using {args.align_threads} thread(s) for Clustal Omega")
    eprint(f"\nRunning Clustal Omega with command:\n {' '.join(clustalo_call)}")
    
    try:
        start_time = time.time()
        subprocess.run(clustalo_call, check=True)
        elapsed = time.time() - start_time
        eprint(f"Processed {input_file} -> {output_file} in {elapsed:.2f}s")
        return output_file # success run
    
    except subprocess.CalledProcessError as e:
        eprint(f"Error during alignment of {input_file}: {e}")
        eprint(f"Command: {' '.join(clustalo_call)}")
        eprint(f"Exit Code: {e.returncode}")
        eprint(f"Error Output: {e.stderr}")
        #raise
        return False  # Return failure flag
    #return output_file


def worker_process(my_file, args, total_workers):
    workerid = int(current_process().name.split("-")[1]) - 1  # Worker ID for debugging
    eprint(f"[Worker {workerid}/{total_workers}] Processing {my_file}")

    try:
        stats = get_sequence_statistics(my_file)
        if args.sequence_stats == "basic":
            eprint(f"Total sequences: {stats['total_sequences']}")
            eprint(f"Minimum length: {stats['min_length']}")
            eprint(f"Maximum length: {stats['max_length']}")
            eprint(f"Average length: {stats['avg_length']:.2f}")
        elif args.sequence_stats == "detailed":
            eprint(f"Total sequences: {stats['total_sequences']}")
            eprint(f"Minimum length: {stats['min_length']}")
            eprint(f"Maximum length: {stats['max_length']}")
            eprint(f"Average length: {stats['avg_length']:.2f}")
            eprint(f"Individual sequence lengths: {stats['sequence_lengths']}")

    except Exception as e:
        eprint(f"[Worker {workerid}/{total_workers}] Failed to read {my_file}: {e}")
        #return None
        return False

    try:
        return get_alignment(my_file, args)
    except Exception as e:
        eprint(f"[Worker {workerid}/{total_workers}] Failed to process {my_file}: {e}")
        return None

    # Call Clustal Omega
    success = get_alignment(my_file, args)
    if not success:
        eprint(f"[Worker {workerid}/{total_workers}] Failed to process {my_file}")
    else:
        eprint(f"[Worker {workerid}/{total_workers}] Successfully processed {my_file}")
    return success    



def main():
    initial_secs = time.time()  # For total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    args = check_args()  # Parse command line arguments

    # Determine input mode and validate
    if args.input_fasta_dir:
        eprint(f" |-- Input mode: Directory ({args.input_fasta_dir})")
        fasta_files = glob.glob(os.path.join(args.input_fasta_dir, f"*{args.extension}"))
        eprint(f" |-- Found {len(fasta_files)} FASTA files in directory.")
    elif args.input_fasta_list:
        eprint(f" |-- Input mode: File List ({args.input_fasta_list})")
        with open(args.input_fasta_list, 'r') as f:
            fasta_files = [line.strip() for line in f if line.strip()]
        eprint(f" |-- Found {len(fasta_files)} FASTA files listed in {args.input_fasta_list}.")
    else:
        exit_with_error("ERROR: No valid input specified. Please provide either --input_fasta_dir or --input_fasta_list.", 1)

    # Track success of all workers
    all_successful = True

    if args.threads > 1:
        with Pool(args.threads) as pool:
            results = list(tqdm(
                pool.starmap(worker_process, [(f, args, args.threads) for f in fasta_files]),
                desc="Processing files in parallel",
                total=len(fasta_files),
            ))
        all_successful = all(results)
    else:
        results = list(tqdm(
            (worker_process(f, args, total_workers=1) for f in fasta_files),
            total=len(fasta_files),
            desc="Processing files sequentially",
        ))
        all_successful = all(results)

    # Only print stats if all files succeeded
    if all_successful:
        output_time(start_time=initial_secs)
        # end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        # total_elapsed_time = elapsed_time(initial_secs)
        # eprint(" '-- Processing complete --'")
        # eprint(f" |-- END TIME: {end_time}")
        # eprint(f" |-- TOTAL ELAPSED TIME: {total_elapsed_time}")
        
    else:
        eprint(" |-- ERROR: One or more files failed to process. Skipping stats.")



    


    
if __name__ == "__main__":
    main()