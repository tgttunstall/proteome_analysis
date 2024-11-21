#!/bin/env bash
#Created on Wed Nov 20 14:51:47 2024
#@author: tanu

if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide as first argument a FASTA file to align"
  >&2 echo "examples:"
#  >&2 echo "   $0 clusters_pig/2.fa"
  >&2 echo "   $0 /path/to/fasta/files/to/align/test.fa"
  exit 22
fi
fasta_file=$1

set -e
set -u

if [ ! -s ${fasta_file} ]; then
  >&2 echo "ERROR: empty or missing ${fasta_file}"
  exit 1
fi

# Extract the base name of the file (without path and extension)
base_name=$(basename ${fasta_file} .fa)

# Define the output file name: output_file would be the shortest match of ${fasta_file}
output_file="${fasta_file%/*}/${base_name}.aln"
sequence_type="Protein"
output_file_format="fa"

# Run Clustal Omega with timing
echo "Starting alignment of ${fasta_file}"
start_time=$(date +%s)

#clustalo -i "$fasta_file" -o "$output_file" --force
#clustalo -i "$fasta_file" -o "$output_file" -t Protein --outfmt fa
#clustalo -i ${fasta_file} -o ${output_file} -t ${sequence_type} --outfmt ${output_file_format}
clustalo -i ${fasta_file} -o ${output_file} -t ${sequence_type} --outfmt ${output_file_format} --force

end_time=$(date +%s)
execution_time=$((end_time - start_time))

# Convert execution time to hours, minutes, and seconds
hours=$((execution_time / 3600))
minutes=$(((execution_time % 3600) / 60))
seconds=$((execution_time % 60))


if [ $? -eq 0 ]; then
  echo "Successfully aligned ${fasta_file}. Output saved to ${output_file}"
  #echo "Alignment completed in ${execution_time} seconds"
  printf "Alignment completed in %02d:%02d:%02d (HH:MM:SS)\n" ${hours} ${minutes} ${seconds}
else
  >&2 echo "ERROR: Clustal Omega alignment failed for $fasta_file"
  exit 1
fi