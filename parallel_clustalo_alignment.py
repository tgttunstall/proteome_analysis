import os
import argparse
from multiprocessing import Pool, current_process
import subprocess
import glob
import time
import sys
from typing import List, Optional
from ffdb import split_file  # Importing split_file from ffdb

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

def get_alignment(input_file, output_dir, clustalo_threads
                  #,force
                  #,sequence_type
                  ):
    """
    Align sequences using Clustal Omega and write to output file.
    """
    base_name = os.path.splitext(os.path.basename(input_file))[0]  # Get base name without extension
    output_file = os.path.join(output_dir, f"{base_name}.aln")  # Construct output file path

    # Call Clustal Omega for alignment
    clustalo_call = [CLUSTALO_EXE,
                     "--infile={}".format(input_file),
                     "--outfmt=fa",
                     "--outfile={}".format(output_file),
                     "--force", 
                     #"--seq_type=Protein",
                     f"--threads={clustalo_threads}"]
    
    #if force:
    #clustalo_call.append("--force")  # Add --force if specified
    
    eprint(f"Running Clustal Omega with command: {' '.join(clustalo_call)}")  # Print the command being run

    start_time = time.time()  # Start timing
    subprocess.run(clustalo_call)
    elapsed = elapsed_time(start_time)  # Calculate elapsed time
    
    eprint(f"Processed {input_file} -> {output_file} in {elapsed}")
    
    return output_file  # Return the name of the output file

def worker_process(my_file, args):
    """Process to work in parallel extracting sequences from FASTA files using Clustal Omega."""
    workerid = int(current_process().name.split("-")[1]) - 1  # 0..threads-1
    if args.progress:
        eprint(f" [{workerid}] working on {my_file}")  # Debug print

    return get_alignment(my_file, args.out_dir, args.clustalo_threads)  # Call get_alignment

def check_args():
    """Parse arguments and check for error conditions."""
    parser = argparse.ArgumentParser(description="Parallel Clustal Omega alignment")
    
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to the input file or directory containing FASTA files.")
    
    parser.add_argument("-o", "--out_dir", type=str, required=True,
                        help="Output directory for aligned files.")
    
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for parallel processing.")
    
    parser.add_argument("-ct", "--clustalo_threads", type=int, default=1,
                        help="Number of threads for Clustal Omega.")
    
    #parser.add_argument("-f", "--force", action="store_true",
    #                    help="Force overwrite of existing output files.")
    
    #parser.add_argument("-st", "--seq_type", type=str, default="Protein",
    #                    help="Sequence type to align.")
    
    parser.add_argument("-e", "--extension", type=str, default=".fa",
                        help="File extension for FASTA files (default: .fa).")
    
    parser.add_argument("--progress", action="store_true",
                        help="Show a progress bar.")

    return parser.parse_args()

def main():
    initial_secs = time.time()  # For total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    
    args = check_args()  # Parse command line arguments
    eprint(f" |-- Input mode: {'Directory' if os.path.isdir(args.input) else 'File List'}")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Determine if input is a directory or a file list
    if os.path.isdir(args.input):
        fasta_files = glob.glob(os.path.join(args.input, f"*{args.extension}"))
        eprint(f" |-- Found {len(fasta_files)} FASTA files in directory.")

        if args.threads > 1:
            with Pool(args.threads) as pool:
                results = pool.starmap(worker_process, [(fasta_file, args) for fasta_file in fasta_files])
            eprint(f" |-- Processed {len(results)} files.")
            eprint(f" |-- Number of workers: {args.threads}")
        else:
            for fasta_file in fasta_files:
                worker_process(fasta_file, args)

    elif os.path.isfile(args.input):
        output_chunk_paths = split_file(args.input)
        eprint(f" |-- Split input into {len(output_chunk_paths)} chunks.")

        if args.threads > 1:
            with Pool(args.threads) as pool:
                results = pool.starmap(worker_process, [(chunk, args) for chunk in output_chunk_paths])
            eprint(f" |-- Processed {len(results)} chunks.")
            eprint(f" |-- Number of workers: {args.threads}")
        else:
            for chunk in output_chunk_paths:
                worker_process(chunk, args)

    else:
        exit_with_error(f"ERROR: The input '{args.input}' is neither a valid directory nor a file.")

    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    total_elapsed_time = elapsed_time(initial_secs)
    eprint(" '-- Processing complete --'")
    eprint(f" |-- END TIME: {end_time}")
    eprint(f" |-- TOTAL ELAPSED TIME: {total_elapsed_time}")

if __name__ == "__main__":
   main()

#    if progress:
#        pbar = tqdm(total=input_file_size, unit="B", position=args.pbar_position)


#python parallel_clustalo_alignment.py -i /home/pub/Work/data_arise_proteome/testC -o /home/pub/Work/data_arise_proteome/testC/results_testC_v0 -t 2 --clustalo_threads 4

python parallel_clustalo_alignment.py -i /home/pub/Work/data_arise_proteome/testC/fasta_list_testC.txt -o /home/pub/Work/data_arise_proteome/testC/results_testC_v0L -t 2 --clustalo_threads 4
