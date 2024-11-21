#!/bin/env bash
#Sun 17 Nov 12:52:18 GMT 2024

if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide as first argument a file containing a list of tsv filepaths"
  >&2 echo "examples:"
  >&2 echo "   $0 ecoli.tsv"
  exit 22
fi
fasta_list=$1

set -e
set -u

if [ ! -s "$fasta_list" ]; then
  >&2 echo "ERROR: empty or missing $fasta_list"
  exit 1
fi

num_jobs=10

SBATCH="${SBATCH/--mail-user=uniprot-prod,$USER/--mail-user=$USER}"
SBATCH="${SBATCH/--mail-type=FAIL,END,/--mail-type=FAIL,} --open-mode=append"

echo "using SBATCH: $SBATCH"

# split the tsv into num_jobs chunks (batches of files)
split -n l/$num_jobs --numeric-suffixes ${fasta_list} ${fasta_list}_CHUNK

pids=""
declare -A childjob
for filebatch in $(ls ${fasta_list}_CHUNK*);
do
  # launch indexer jobs to index all files in the batch
  $SBATCH --wait -J "indexing fasta files in $filebatch" -o ${filebatch}.out --wrap="./fasta_indexer.sh $filebatch" &
  pids="${!} ${pids}"
  childjob+=([${!}]="fasta indexer for files in ${filebatch}; check ${filebatch}.out")
done
jobsfailed=0
for pid in $pids
do
  wait $pid
  return_code=$?
  if [ "${return_code}" -ne 0 ]; then
    >&2 echo "    =>ERROR: exit code ${return_code} for ${childjob[${pid}]}"
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
