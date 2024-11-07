#!/usr/bin/env python3
# changelog
# Wed  6 Nov 2024 10:05:27 GMT 0.1 started
# Wed  6 Nov 2024 10:05:27 GMT 1.0 working sequential approach

# imports
import os
import sys
import re
import time
import argparse
from glob import glob

# constants
CHUNKSIZE = 1024 * 1024  # 1Mb
RE_REMOVE_EXTENSION = re.compile(r"\.[^.]+$")
DESCRIPTION = """
Script to label the second column of a tsv file,
assumed to contain protein id from fasta headers,
with the name of the fasta file where the proteins are located.
It will also compact the first column of the tsv file,
assumed to be cluster identifiers, into a list of sequential integer numbers.
An additional column will be added to mark the first protein_id of each cluster
(assumed to be the representative protein of the cluster)
Order of the output file will be the same as the input file

Arguments needed:
    FASTA_DIR Directory containing the fasta files whose headers we will match
    TSV_FILE The input file (second column assumed to have protein identifiers)
    OUT_FILE The output file that will be created

Sample input file:
    ENSSSCP00000055324|661\tENSSSCP00000055324|661
    ENSSSCP00000055324|661\tENSSSCP00055011301|568
    ENSSSCP00000055324|661\tENSSSCP00035021320|596

Sample output file:
    cluster_id\tprotein_id\tproteomes\tis_rep
    0\tENSSSCP00000055324|661\tproteome_35497\t*
    0\tENSSSCP00055011301|568\tproteome_4698922\t
    0\tENSSSCP00035021320|596\tproteome_4698918\t

Example call:
    ./label_clusters.py --fasta_dir pig/ --tsv_file <(uniq results_pig/Specie_protein_cluster.tsv) --out_file results_pig/Labelled_Specie_protein_cluster.tsv --prefix proteome_ --extension .fa
"""


# helper functions
def secs2time(secs):
    """
    human readable printout of seconds elapsed
    e.g.:
    secs2time(3663)
    """
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)


def elapsed_time(start_time, work_done=None):
    """
    compute elapsed time from given start_time (in seconds) and return
    well formatted string. If work_done argument given, compute speed of the process
    e.g.:
    start_secs = time.time(); iterations_done = 10; time.sleep(2); print(" '-- Elapsed: {}, {} it/s --'".format(*elapsed_time(start_secs, iterations_done)))
    print(" '-- Elapsed: {} --'".format(elapsed_time(start_secs)))
    """
    process_time = time.time() - start_time
    if work_done is None:
        return secs2time(process_time)
    process_speed = round(work_done / process_time, 2)
    return secs2time(process_time), process_speed


def eprint(*myargs, **kwargs):
    """
    print to stderr, useful for error messages and to not clobber stdout
    """
    print(*myargs, file=sys.stderr, **kwargs)


def read_file_in_chunks(file_path, chunk_size=CHUNKSIZE):
    """
    process a file in chunks yielding only full lines
    """
    with open(file_path, "r") as file:
        buffer = ""
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
                break  # end of file

            buffer += chunk
            lines = buffer.splitlines(True)

            # process all lines except the last one which may be incomplete
            for line in lines[:-1]:
                yield line.strip()

            # keep the remaining possibly incomplete line in the buffer
            buffer = lines[-1]

        if buffer:
            yield buffer.strip()


def delete_files(filenames, path=None):
    """
    delete temporary sub-files
    optionally specify a path where to look for the files
    """
    for filename in filenames:
        if path is not None:
            filename = os.path.join(path, filename)
        if os.path.isfile(filename):
            os.remove(filename)


def print_stats(start_time, protein_id_count, cluster_count):
    """
    print some final statistics on number of identifiers mapped and time
    """
    eprint(f" '-- {protein_id_count} protein identifiers with {cluster_count} clusters")
    eprint(
        " '-- Elapsed: {}, {} ids/sec --'".format(
            *elapsed_time(start_time, protein_id_count)
        )
    )


# functions
def check_args(DESCRIPTION):
    """
    parse arguments and check for error conditions
    """

    class CustomArgumentParser(argparse.ArgumentParser):
        def print_help(self, *args, **kwargs):
            """
            print custom text before the default help message
            """
            print(DESCRIPTION)
            super().print_help(*args, **kwargs)

    parser = CustomArgumentParser(
        description="Proteome labeller for lists of protein indentifiers."
    )
    parser.add_argument(
        "--fasta_dir",
        type=str,
        required=True,
        help="Directory containing all .fa files",
    )
    parser.add_argument(
        "--tsv_file",
        type=argparse.FileType("r"),
        required=True,
        help="Path to the input file",
    )
    parser.add_argument(
        "--out_file", type=str, required=True, help="Path to the output file"
    )
    parser.add_argument(
        "--extension",
        type=str,
        required=False,
        help="Extension for files in fasta_dir, e.g. '.fa'",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        required=False,
        help="Optionally prefix for filenames, to be removed; e.g. 'proteome_'",
    )
    parser.add_argument(
        "--nolabel",
        dest='unmapped_label',
        type=str,
        required=False,
        help="Optional character or string to use when no mapping could be found; e.g. '?'",
        default='',
    )
    args = parser.parse_args()
    # check for out_file presence and writing ability
    if os.path.isfile(args.out_file):
        eprint(f"    => ERROR: file '{args.out_file}' exists")
        eprint("       please remove it as we refuse to overwrite it!")
        sys.exit(1)
    try:
        myoutputfh = open(args.out_file, "w")
        myoutputfh.close()
    except PermissionError:
        eprint(f"    => ERROR: Cannot open file '{args.out_file}' for writing")
        sys.exit(1)  # permission
    delete_files([args.out_file])  # delete the just touched file

    # check for fasta_dir
    if args.extension is not None:  # get only files with specified extension
        if args.prefix is not None:
            fasta_files = glob(
                os.path.join(args.fasta_dir, f"{args.prefix}*{args.extension}")
            )
        else:
            fasta_files = glob(os.path.join(args.fasta_dir, f"*{args.extension}"))
        if not fasta_files:
            eprint(
                f"No '{args.extension}' files found in {args.fasta_dir}. Please check the directory."
            )
            sys.exit(2)  # no such file or dir
    else:  # get all files
        fasta_files = glob(os.path.join(args.fasta_dir, f"*"))
        if not fasta_files:
            eprint(f"No files found in {args.fasta_dir}. Please check the directory.")
            sys.exit(2)  # no such file or dir

    return args, fasta_files


def create_proteome_protein_map(fasta_files, read_method="lines"):
    """
    process all fasta_files sequentially to create a dictionary protein_id: proteome
    three different reading methods are provided: 'chunks' 'lines' 'full'
    performance considerations:
        it takes 1s to process 13 files of approx 30M each for a total of 627957 ids
        'chunks': ~800000 ids/s
        'full': ~770000 ids/s
        'lines': ~1100000 ids/s
        (tested on a linux datacentre)
        (on macos 'full' is as fast as 'lines' and 'chunks' slightly less)
    'chunks' may or may not be useful for very big fasta files, to be tested

    """
    proteome_protein_map = {}

    def _process_protein_id(line):
        """
        store protein_id extracted from fasta headers
        """
        protein_id = line.rstrip()[1:]
        if protein_id not in proteome_protein_map:
            proteome_protein_map[protein_id] = proteome_id
        else:
            proteome_protein_map[protein_id] += " " + proteome_id

    for file in fasta_files:
        proteome_id = os.path.basename(file)

        # remove extension, if any
        if args.extension is not None:
            proteome_id = proteome_id[0 : -len(args.extension)]
        else:
            proteome_id = RE_REMOVE_EXTENSION.sub("", proteome_id)

        # optionally remove prefix from the filenames
        if args.prefix is not None:
            proteome_id = proteome_id[len(args.prefix) :]

        if read_method == "lines":  # read line by line
            with open(file, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        _process_protein_id(line)
        elif read_method == "full":  # read whole file in memory
            with open(file, "r") as f:
                for line in f.readlines():
                    if line.startswith(">"):
                        _process_protein_id(line)
        elif read_method == "chunks":  # read file in chunks
            for line in read_file_in_chunks(file):
                if line.startswith(">"):  # fasta header
                    _process_protein_id(line)
        else:
            eprint(f"    => ERROR: no such read_method {read_method}")

    # eprint(list(proteome_protein_map.items())[:5]) #debug, first 5 items in the map
    # eprint(list(proteome_protein_map.items())[-5:]) #debug, last 5 items in the map
    return proteome_protein_map


def label_proteins(input_fh, output_file, proteome_protein_map, unmapped_label=''):
    """
    use the provided mapping to label the input_file, assigning proteome labels to protein_ids
    """
    cluster_counter = -1
    prev_cluster = None

    with open(output_file, "w") as outfile:
        outfile.write("cluster_id\tprotein_id\tproteomes\tis_rep\n")
        for line in input_fh:
            cluster_id, protein_id = line.rstrip().split("\t")
            if cluster_id != prev_cluster:
                prev_cluster = cluster_id
                cluster_counter += 1
                outfile.write(
                    f"{cluster_counter}\t{protein_id}\t{proteome_protein_map.get(protein_id, unmapped_label)}\t*\n"
                )
            else:
                outfile.write(
                    f"{cluster_counter}\t{protein_id}\t{proteome_protein_map.get(protein_id, unmapped_label)}\t\n"
                )

    return cluster_counter + 1


if __name__ == "__main__":
    # ===============
    # 0: argument parsing and input/output files checking
    # ===============
    initial_secs = time.time()  # for total time count
    args, fasta_files = check_args(DESCRIPTION)
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" .-- BEGUN {timestamp} --. ")

    # ===============
    # 1: create protein->proteome map
    # ===============
    start_secs = time.time()
    proteome_protein_map = create_proteome_protein_map(fasta_files, read_method="lines")
    eprint(
        f" |-- processed {len(proteome_protein_map)} protein ids from {len(fasta_files)} proteome files"
    )
    eprint(
        " |-- mapping completed {} -- Elapsed: {}, {} ids/s --".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            *elapsed_time(start_secs, len(proteome_protein_map)),
        )
    )

    # ===============
    # 2: assign proteome labels to protein identifiers
    # ===============
    start_secs = time.time()
    clusters_count = label_proteins(args.tsv_file, args.out_file, proteome_protein_map, args.unmapped_label)
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" |-- processed {clusters_count} clusters")
    eprint(
        " |-- labelling completed -- Elapsed: {}, {} ids/s --".format(
            *elapsed_time(start_secs, len(proteome_protein_map)),
        )
    )

    # total time
    eprint(
        " |-- total time: {} --".format(
            elapsed_time(initial_secs),
        )
    )
    eprint(f" .-- ENDED {timestamp} --. ")
