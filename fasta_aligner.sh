#!/bin/env bash
#Created on Wed Nov 20 17:28:50 2024
#@author: tanu

if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide a directory containing FASTA files to align"
  >&2 echo "example:"
  >&2 echo "   $0 /path/to/fasta/dir"
  exit 22
fi
fasta_dir=$1

set -e
set -u

if [ ! -d "${fasta_dir}" ]; then
  >&2 echo "ERROR: Directory ${fasta_dir} does not exist or is not a directory"
  exit 1
fi

echo "FASTA directory: ${fasta_dir}"

# Count the number of FASTA files
fasta_count=$(find "${fasta_dir}" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) | wc -l)
echo "Number of FASTA files in ${fasta_dir}: ${fasta_count}"

sequence_type="Protein"
output_file_format="fa"
n_threads=8

printf "\nStarting...\n\n"
start_time_whole=$(date +%s)

current_file=0
for fasta_file in ${fasta_dir}/*.fa ${fasta_dir}/*.fasta; do
  if [ ! -f ${fasta_file} ]; then
    continue
  fi
  
  # Count the number of sequences using grep
  seq_count=$(grep -c "^>" "${fasta_file}")
  #seq_count=$(awk 'BEGIN {count=0} /^>/ {count++} END {print count}' "${fasta_file}")

  ((current_file+=1))
  #printf "\nProcessing ${current_file}/${fasta_count}: ${fasta_file}"
  printf "\nProcessing ${current_file}/${fasta_count}"
  
  #printf "\nStarting alignment of ${fasta_file}, \nNumber of sequences = ${seq_count}"
  printf "\nStarting alignment: ${fasta_file}"

  base_name=$(basename ${fasta_file} .fa)
  base_name=${base_name%.fasta}
  output_file="${fasta_dir}/${base_name}.aln"

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
printf "\n"${fasta_count}" fasta file alignments completed in %02d:%02d:%02d (HH:MM:SS)\n" ${hours_whole} ${minutes_whole} ${seconds_whole}

##########
#foo=0

#until [ $foo -gt 99 ]
#do
#  echo foo: $foo
##  ((foo+=1))
#  ((foo=foo+1))
#done