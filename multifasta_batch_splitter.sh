#!/bin/bash
#input_fasta="/home/pub/Work/data_arise_proteome/ggcaller/gene_calls.faa""
#num_files=2  # Number of output files to create
#
## Count the number of sequences (header lines starting with ">")
#total_seqs=$(grep -c "^>" $input_fasta)
#seqs_per_file=$(( (total_seqs + num_files - 1) / num_files )) # Ceiling division
#
## Split the file
#awk -v n=$seqs_per_file -v prefix="split_" '
#  /^>/ { if (count % n == 0) file = prefix ++file_idx ".fasta"; count++ }
#  { print >> file }
#' $input_fasta

# Usage function
usage() {
    echo "Usage: $0 -i <input_fasta> -n <num_files> -o <output_dir>"
    exit 1
}

# Parse command-line arguments
while getopts "i:n:o:" opt; do
    case $opt in
        i) input_fasta="$OPTARG" ;;
        n) num_files="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure required arguments are provided
if [[ -z "$input_fasta" || -z "$num_files" || -z "$output_dir" ]]; then
    usage
fi

# Validate input file
if [[ ! -f "$input_fasta" ]]; then
    echo "Error: Input file $input_fasta does not exist."
    exit 1
fi

# Validate and create output directory
mkdir -p "$output_dir" || { echo "Error: Cannot create output directory $output_dir."; exit 1; }

# Count the number of sequences (header lines starting with ">")
total_seqs=$(grep -c "^>" "$input_fasta")
seqs_per_file=$(( (total_seqs + num_files - 1) / num_files )) # Ceiling division

# Split the file
awk -v n="$seqs_per_file" -v prefix="$output_dir/split_" '
    /^>/ {
        if (count % n == 0) file = prefix file_idx++ ".fasta"; 
        count++;
    }
    { print >> file }
' "$input_fasta"

echo "Split completed. Files are saved in $output_dir"
