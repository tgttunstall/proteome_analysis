#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:16:53 2024

@author: tanu
"""
import os
import argparse
from multiprocessing import Pool, current_process, cpu_count
import subprocess
import glob
import time
import sys
from typing import List, Optional
from tqdm import tqdm

################
BUFFERSIZE = 1048576  # 1Mb
DESCRIPTION = """
Script to parallelise clustal omega alignments
[a little more documentation, example calls, etc]
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


def elapsed_time(start_time, work_done=None, invert=False):
    """
    Computes the elapsed time from a given start time in seconds and returns a formatted string.
    If `work_done` is specified, also computes the speed of the process.

    Args:
        start_time (float): The start time in seconds (from time.time()).
        work_done (int, optional): Number of completed iterations or tasks.

    Returns:
        str or tuple: A formatted string with elapsed time if `work_done` is None,
                      otherwise a tuple with formatted elapsed time and computed speed in "it/s".

    Example:
        start_secs = time.time()
        time.sleep(2)
        print(" '-- Elapsed: {} --'".format(elapsed_time(start_secs)))
        # Output example: " '-- Elapsed: 00h 00m 02s --'"

        iterations_done = 10
        print(" '-- Elapsed: {}, {} it/s --'".format(*elapsed_time(start_secs, iterations_done)))
        # Output example: " '-- Elapsed: 00h 00m 02s, 5.0 it/s --'"
    """
    process_time = time.time() - start_time
    if work_done is None:
        return secs2time(process_time)
    if invert:
        process_speed = round(process_time / work_done, 2)
    else:
        process_speed = round(work_done / process_time, 2)
    return secs2time(process_time), process_speed


def exit_with_error(message: str, code: int = 1):
    """Prints an error message to stderr and exits the program with the specified exit code."""
    eprint(f" => {message}")
    sys.exit(code)


def eprint(*myargs, **kwargs):
    """Prints the provided arguments to stderr."""
    print(*myargs, file=sys.stderr, **kwargs)


def delete_files(filenames: List[str], path: Optional[str] = None):
    """
    Deletes specified temporary files.

    Args:
        filenames (List[str]): List of filenames to delete.
        path (Optional[str]): Optional base directory to prepend to each filename.

    Returns:
        None
    """
    for filename in filenames:
        if path is not None:
            filename = os.path.join(path, filename)
        if os.path.isfile(filename):
            os.remove(filename)

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
    exclusive_group = input_group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument("-d", "--fasta_dir", type=str, required=False,
                                 dest="input_fasta_dir",
                                 help="Directory containing fasta files to be aligned.")
    exclusive_group.add_argument("-l", "--fasta_list", type=is_valid_file, required=False,
                                 dest="input_fasta_list",
                                 help="Path to a text file containing a list of FASTA files.")

    # General options
    parser.add_argument("-t", "--threads", type=positive_integer, default=1,
                        help="Number of threads for parallel processing.")
    parser.add_argument("-e", "--extension", type=str, default=".fa",
                        help="File extension for FASTA files (default: .fa).")

    parser.add_argument("--sequence_stats", choices=["none", "basic", "detailed"], default="none",
                        help="Level of sequence statistics to report (none, basic, or detailed)")

    #==================================
    # Clustal Omega-specific arguments
    #==================================

    clustalo_group = parser.add_argument_group("Clustal Omega specific arguments")

    clustalo_group.add_argument("-o", "--out_dir", type=str, required=True,
                                help="Output directory for aligned files.")

    clustalo_group.add_argument("-at", "--align_threads", type=positive_integer, default=1,
                                help="Number of threads for Clustal Omega.")

    clustalo_group.add_argument("-f", "--force", action="store_true",
                                help="Force overwrite of existing output files.")

    clustalo_group.add_argument("-st", "--seqtype", type=str, choices=["Protein", "DNA"],
                                default="Protein",
                                help="Specify the sequence type (e.g., 'Protein' or 'DNA').")

    clustalo_group.add_argument("-of", "--outfmt", type=str,
                                choices=["fa","clu","msf","phy","selex","st","vie"], default="fa",
                                help="Specify the MSA output format.")

    args = parser.parse_args()

    # Validate input arguments
    #if args.input_fasta_dir and args.input_fasta_list:
    #   exit_with_error("ERROR: Please specify either --input_fasta_dir or --input_fasta_list, not both.", 1)

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

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")

    return args


def get_sequence_statistics(fasta_file):
    """[docstring here]"""
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

    if not os.path.isfile(input_file):
        eprint(f" |-- ERROR: no such file {input_file}")
        return 0

    # Construct Clustal Omega command
    clustalo_call = [
        CLUSTALO_EXE,
        "--infile", input_file,   #f"--infile={input_file}",
        "--outfile", output_file, #"--outfile", "{}".format(output_file)
        "--outfmt", args.outfmt,
        "--threads", str(args.align_threads), #because subprocess.run expects it!
        "--seqtype", args.seqtype,
    ]

    if os.path.isfile(output_file):
        if args.force:
            clustalo_call.append("--force")
        else:
            eprint(f" |-- WARNING: Cowardly refusing to overwrite already existing file '{output_file}'. Use --force to force overwriting.")
            return 0

    # Log the command and the number of threads being used
    #eprint(f"\nForce overwrite enabled: {args.force}")

    #eprint(f"Using {args.align_threads} thread(s) for Clustal Omega") #debug
    #eprint(f"\nRunning Clustal Omega with command:\n {' '.join(clustalo_call)}") #debug

    try:
        #start_time = time.time()
        subprocess.run(clustalo_call, check=True)
        #elapsed = time.time() - start_time
        #eprint(f"Processed {input_file} -> {output_file} in {elapsed:.2f}s") #debug
        return 1 # success run

    except subprocess.CalledProcessError as e:
        eprint(f"Error during alignment of {input_file}")
        eprint(f"Command: {' '.join(clustalo_call)}")
        eprint(f"Exit Code: {e.returncode}")
        eprint(f"Error Output: {e.stderr}")
        #raise
        return 0 # Return failure flag
    #return output_file


def initializer(args_arg):
    """
    initializer to set global variables for workers
    """
    global args
    args = args_arg


def worker_process(my_file):
    """
    Process a single file for alignment and optionally write sequence statistics.
    TODO: Add more
    """
    workerid = int(current_process().name.split("-")[1]) - 1  # Worker ID for debugging
    #eprint(f"[Worker {workerid}/{args.threads}] Processing {my_file}") #debug
    
    tmpstat_filename = f"tmpstat{workerid}"
    tmpstat_content = []
    
    if args.sequence_stats != "none":
        try:
            stats = get_sequence_statistics(my_file)
            #with open(f"tmpstat{workerid}", 'a') as statfh:
                #statfh.write(....)
            if args.sequence_stats == "basic":
                tmpstat_content.append(f"File: {my_file}\n")
                tmpstat_content.append(f"Total sequences: {stats['total_sequences']}\n")
                tmpstat_content.append(f"Minimum length: {stats['min_length']}\n")
                tmpstat_content.append(f"Maximum length: {stats['max_length']}\n")
                tmpstat_content.append(f"Average length: {stats['avg_length']:.2f}\n")
            elif args.sequence_stats == "detailed":
                tmpstat_content.append(f"File: {my_file}\n")
                tmpstat_content.append(f"Total sequences: {stats['total_sequences']}\n")
                tmpstat_content.append(f"Minimum length: {stats['min_length']}\n")
                tmpstat_content.append(f"Maximum length: {stats['max_length']}\n")
                tmpstat_content.append(f"Average length: {stats['avg_length']:.2f}\n")
                tmpstat_content.append(f"Individual sequence lengths: {stats['sequence_lengths']}\n")

            with open(tmpstat_filename, 'a') as statfh:
                statfh.writelines(tmpstat_content)

        except Exception as e:
            eprint(f"[Worker {workerid}/{args.threads}] Failed to read {my_file}: {e}")
            return False

    # Directly return the result of get_alignment
    return get_alignment(my_file, args)    

if __name__ == "__main__":
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

    fasta_files_count = len(fasta_files)
    results = []
    
    tempstat_files = [f"tmpstat{workerid}" for workerid in range(args.threads)]
    with Pool(processes=args.threads, initializer=initializer, initargs=(args,)) as pool:
        for result in tqdm(
            pool.imap_unordered(worker_process, fasta_files),
            desc="aligning",
            total=fasta_files_count,
        ):
            results.append(result)
    successful_count = sum(results)

    #if args.stats....
    #tempstatfiles = [tmpstat + workerid for workerid in range(args.threads)]
    #with open(final_stat_file, 'w') as finalstatsfh:
    #     for tempstatfile in tempstatfiles:
    #     with open(tempstatfile, 'r') as statfh:
    #         finalstatsfh.write(statfh.read())
    #delete_files(tempstatfiles)
    
    if args.sequence_stats != "none":
        final_stats_file = os.path.join(args.out_dir, "final_stats.txt")
        with open(final_stats_file, 'w') as finalstatsfh:
            for tempstat_file in tempstat_files:
                if os.path.exists(tempstat_file):
                    with open(tempstat_file, 'r') as statfh:
                        finalstatsfh.write(statfh.read())
        delete_files(tempstat_files)
        eprint(f"Final sequence stats written to {final_stat_file}")
        
    # Only print stats if all files succeeded
    if successful_count != fasta_files_count:
        eprint(f" |-- ERROR: {fasta_files_count - successful_count} file(s) failed to process.")

    if successful_count == 0:
        eprint(" |-- ERROR: no alignment created")
    else:
        eprint("|-- {} alignments completed -- Elapsed: {}, {} s/alignment --".format(successful_count, *elapsed_time(initial_secs + 1E-10, successful_count, invert=True),))
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
