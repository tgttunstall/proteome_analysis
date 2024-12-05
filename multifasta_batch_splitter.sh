#!/bin/bash
# Usage function
usage() {
    echo "Usage: $0 -i <input_fasta> -n <num_files> -o <output_dir> [-p]"
    echo "  -i: Input multi-FASTA file"
    echo "  -n: Number of output files"
    echo "  -o: Output directory"
    echo "  -p: Enable progress reporting (optional)"
    exit 1
}

# Parse command-line arguments
progress=false
while getopts "i:n:o:p" opt; do
    case $opt in
        i) input_fasta="$OPTARG" ;;
        n) num_files="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        p) progress=true ;;
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
if [[ $total_seqs -eq 0 ]]; then
    echo "Error: No sequences found in $input_fasta."
    exit 1
fi
echo "Total sequences found: $total_seqs"

seqs_per_file=$(( (total_seqs + num_files - 1) / num_files )) # Ceiling division
echo "Splitting into $num_files files with approximately $seqs_per_file sequences per file."

outfname=$(basename ${input_fasta})                                                                 
outfname=${outfname}_split_                                                                         
                                                                                                     
# Split the file                                                                                    
#awk -v n="$seqs_per_file" -v prefix="$output_dir/split_" -v total="$total_seqs" -v progress="$progress" '
awk -v n="$seqs_per_file" -v prefix="${output_dir}/${outfname}" -v total="$total_seqs" -v progress="$progress" '
    BEGIN { file_idx = 1; count = 0; print "Starting split..."; }
    /^>/ { 
        if (count % n == 0) {
            if (count > 0) close(file);
            file = prefix file_idx++ ".fasta"; 
            printf "Creating file: %s\n", file > "/dev/stderr";
        }
        count++;
        if (progress == "true" && (count % 100 == 0 || count == total)) {
            printf "Progress: %d/%d sequences processed...\n", count, total > "/dev/stderr";
        }
    }
    { print >> file }
    END {
        printf "Finished splitting into %d files.\n", file_idx - 1 > "/dev/stderr";
    }
' "$input_fasta"

echo "Split completed. Files are saved in $output_dir"
