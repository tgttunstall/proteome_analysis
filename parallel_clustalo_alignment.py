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
        "--threads", str(args.align_threads), #because subprocess.run excepts it!
        "--seqtype", args.seqtype,
    ]

    # Add --force if specified
    if args.force:
        clustalo_call.append("--force")

    # Log the command and the number of threads being used
    eprint(f"Using {args.align_threads} thread(s) for Clustal Omega")
    #eprint(f"\nForce overwrite enabled: {args.force}")
    eprint(f"\nRunning Clustal Omega with command:\n {' '.join(clustalo_call)}")

    start_time = time.time()

    try:
        subprocess.run(clustalo_call, check=True)
        elapsed = time.time() - start_time
        eprint(f"Processed {input_file} -> {output_file} in {elapsed:.2f}s")
    except subprocess.CalledProcessError as e:
        eprint(f"Error during alignment of {input_file}: {e}")
        raise

    return output_file

def worker_process(my_file, args, total_workers):
    """
    Process a single file in parallel for Clustal Omega.
    """
    workerid = int(current_process().name.split("-")[1]) - 1  # Worker ID for debugging
    eprint(f"[Worker {workerid}/{total_workers}] Processing {my_file}")

    try:
        return get_alignment(my_file, args)
    except Exception as e:
        eprint(f"[Worker {workerid}/{total_workers}] Failed to process {my_file}: {e}")
        return None

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

        if args.threads > 1:
            with Pool(args.threads) as pool:
                results = list(tqdm(
                    pool.starmap(worker_process, [(fasta_file, args, args.threads) for fasta_file in fasta_files]),
                    desc="Processing files in parallel",
                    total=len(fasta_files)
                ))
            eprint(f" |-- Processed {len(results)} files.")
        else:
            for fasta_file in tqdm(fasta_files, desc="Processing files sequentially"):
                worker_process(fasta_file, args, total_workers=1)

    elif args.input_fasta_list:
        eprint(f" |-- Input mode: File List ({args.input_fasta_list})")
        
        # Read the list of FASTA files from the provided file
        with open(args.input_fasta_list, 'r') as f:
            fasta_files = [line.strip() for line in f if line.strip()]  # Read lines and strip whitespace

        eprint(f" |-- Found {len(fasta_files)} FASTA files listed in {args.input_fasta_list}.")

        if args.threads > 1:
            with Pool(args.threads) as pool:
                results = list(tqdm(
                    pool.starmap(worker_process, [(fasta_file, args, args.threads) for fasta_file in fasta_files]),
                    desc="Processing files in parallel",
                    total=len(fasta_files)
                ))
            eprint(f" |-- Processed {len(results)} files.")
        else:
            for fasta_file in tqdm(fasta_files, desc="Processing files sequentially"):
                worker_process(fasta_file, args, total_workers=1)

    else:
        exit_with_error("ERROR: No valid input specified. Please provide either --input_fasta_dir or --input_fasta_list.", 1)

    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    total_elapsed_time = elapsed_time(initial_secs)
    eprint(" '-- Processing complete --'")
    eprint(f" |-- END TIME: {end_time}")
    eprint(f" |-- TOTAL ELAPSED TIME: {total_elapsed_time}")

    
if __name__ == "__main__":
    main()

###############################################################################
#python parallel_clustalo_alignment.py -d /home/pub/Work/data_arise_proteome/testC -o /home/pub/Work/data_arise_proteome/testC/results_testC_v0L -t 2 -at 4 --force
#python parallel_clustalo_alignment.py -f /home/pub/Work/data_arise_proteome/testC/fasta_list_testC.txt -o /home/pub/Work/data_arise_proteome/testC/results_testC_v0L -t 2 -at 4 --force
