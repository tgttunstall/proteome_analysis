#!/bin/env bash
#Created on Wed Nov 20 17:28:50 2024
#@author: tanu

if [ $# -lt 2 ]; then
  >&2 echo "Usage: $0 <input_fasta_dir> <output_dir> [threads]"
  >&2 echo "Example:"
  >&2 echo "   $0 /path/to/fasta/dir /path/to/output/dir 8"
  exit 22
fi

fasta_dir=$1
output_dir=$2
n_threads=${3:-8}  # Default to 8 threads if not specified

set -e
set -u

if [ ! -d "${fasta_dir}" ]; then
  >&2 echo "ERROR: Input directory ${fasta_dir} does not exist or is not a directory"
  exit 1
fi

if [ ! -d "${output_dir}" ]; then
  >&2 echo "Creating output directory: ${output_dir}"
  mkdir -p "${output_dir}"
fi

echo "FASTA directory: ${fasta_dir}"
echo "Output directory: ${output_dir}"
echo "Number of threads: ${n_threads}"

# Count the number of FASTA files
fasta_count=$(find "${fasta_dir}" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) | wc -l)
echo "Number of FASTA files in ${fasta_dir}: ${fasta_count}"

sequence_type="Protein"
output_file_format="fa"

printf "\nStarting...\n\n"
start_time_whole=$(date +%s)

current_file=0
for fasta_file in ${fasta_dir}/*.fa ${fasta_dir}/*.fasta; do
  if [ ! -f ${fasta_file} ]; then
    continue
  fi
  
  seq_count=$(grep -c "^>" "${fasta_file}")

  ((current_file+=1))
  printf "\nProcessing ${current_file}/${fasta_count}"
  printf "\nStarting alignment: ${fasta_file}"

  base_name=$(basename ${fasta_file} .fa)
  base_name=${base_name%.fasta}
  output_file="${output_dir}/${base_name}.aln"

  printf "\nNumber of sequences = ${seq_count}"
  start_time=$(date +%s)

  clustalo -i "${fasta_file}" -o "${output_file}" -t "${sequence_type}" --outfmt "${output_file_format}" --threads ${n_threads} --force

  end_time=$(date +%s)
  execution_time=$((end_time - start_time))

  hours=$((execution_time / 3600))
  minutes=$(((execution_time % 3600) / 60))
  seconds=$((execution_time % 60))

  if [ $? -eq 0 ]; then
    printf "\nSuccessfully aligned: ${fasta_file}. \nOutput saved: ${output_file}\n"
    printf "\nAlignment completed in %02d:%02d:%02d (HH:MM:SS)\n" $hours $minutes $seconds
  else
    >&2 printf "\nERROR: Clustal Omega alignment failed for ${fasta_file}"
  fi
  
  echo "----------------------------------------"

done
end_time_whole=$(date +%s)
execution_time_whole=$((end_time_whole - start_time_whole))

hours_whole=$((execution_time_whole / 3600))
minutes_whole=$(((execution_time_whole % 3600) / 60))
seconds_whole=$((execution_time_whole % 60))

printf "\nAll alignments completed."
printf "\n${fasta_count} fasta file alignments completed in %02d:%02d:%02d (HH:MM:SS)\n" ${hours_whole} ${minutes_whole} ${seconds_whole}
