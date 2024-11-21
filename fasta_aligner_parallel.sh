#!/bin/bash
#Created on Thu Nov 21 18:14:39 2024
#@author: tanu


if [ "X$1" == "X" ] || [ "X$2" == "X" ]; then
  >&2 echo "Usage: $0 <fasta_list_file> <output_directory>"
  >&2 echo "Example:"
  >&2 echo "   $0 fasta_files.txt /path/to/output/dir"
  exit 22
fi

fasta_list=$1
output_dir=$2

set -e
set -u

if [ ! -s "$fasta_list" ]; then
  >&2 echo "ERROR: empty or missing $fasta_list"
  exit 1
fi

if [ ! -d "$output_dir" ]; then
  >&2 echo "Creating output directory: $output_dir"
  mkdir -p "$output_dir"
fi

num_jobs=10
threads_per_job=15

SBATCH="sbatch --parsable --mail-user=$USER --mail-type=FAIL --open-mode=append"

echo "Using SBATCH: $SBATCH"

# Split the fasta list into num_jobs chunks (batches of files)
split -n l/$num_jobs --numeric-suffixes ${fasta_list} ${fasta_list}_CHUNK

pids=""
declare -A childjob

for filebatch in ${fasta_list}_CHUNK*; do
  # Launch aligner jobs to process all files in the batch
     job_id=$($SBATCH --cpus-per-task=$threads_per_job \
           --job-name="align_batch_$(basename $filebatch)" \
           --output="${filebatch}.out" \
           --wrap="./fasta_aligner.sh $filebatch $output_dir $threads_per_job")
  
  pids="$job_id $pids"
  childjob+=([${job_id}]="fasta aligner for files in ${filebatch}; check ${filebatch}.out")
done

echo "Waiting for all jobs to complete..."
jobsfailed=0
for job_id in $pids; do
  scontrol wait job $job_id
  return_code=$?
  if [ "${return_code}" -ne 0 ]; then
    >&2 echo "    =>ERROR: exit code ${return_code} for ${childjob[${job_id}]}"
    jobsfailed=$((jobsfailed + 1))
  fi
done

if (( jobsfailed == 0 )); then
  echo "ALL DONE, cleaning up intermediate files and logs"
  rm ${fasta_list}_CHUNK*.out
  rm ${fasta_list}_CHUNK*[0-9]
else
  >&2 echo "    =>$jobsfailed job(s) failed, you may want to redo them singularly or restart the whole process"
fi
