#!/bin/env bash
#Wed  9 Oct 17:45:48 BST 2024 1.0
#Tue 26 Nov 15:29:16 GMT 2024 1.1 added parameter for covmode
#Mon  2 Dec 13:06:12 GMT 2024 1.2 adapted for use with standard mmseqs

#NOTE: this is wrapper for the standard mmseqs (from https://github.com/soedinglab/MMseqs2)

set -e
set -u

#CONSTANTS and SETTINGS
#paths:
BIGFASTAPATH=$WUP/temp/mmseqs/fasta #path where to store concatenated proteome fasta files
MMSEQS_EXEC="./mmseqs" #path to mmseqs executable
MMSEQS_EXEC="/hps/software/users/jlees/tanushree/MMseqs2/build/bin/mmseqs"

#parameters for clustering:
seq_id_mode=0 #--seq-id-mode 0: alignment length 1: shorter, 2: longer sequence
#seq-id-mode 0: calculate sequence identity as number of identical residues divided by the length of the aligned region
#seq-id-mode 1: calculate sequence identity as identical AA match / min(qLen, tLen)

#parameters for resources:
threads=32 #number of parallel threads
mem="500G" #memory
time="7-0" #requested execution time

#ARGUMENT PARSING
extension=".fa" #default
fasta_path=""
coverage=""
min_seq_id=""
cov_mode=""
usage="Usage: $0 [-hd:e:m:c:p:]
   -h   Show this help message
   -d   Directory of (proteome) fasta files to cluster
   -e   Extension for fasta files; default: '.fa'
   -m   Covmode (--cov-mode): 0) bidirectional, 1) target coverage, 2) query coverage
   -c   Minimum coverage (-c): [0.0-1.0] as the # of aligned residue pairs divided by... (according to covmode)
   -p   Minimum sequence identity (--min-seq-id): [0.0-1.0] as the # of identical aligned residues divided by the number of aligned columns including internal gap columns"

while getopts ":hd:e:m:c:p:" opt; do
  case $opt in
    h)
      echo "$usage"
      exit 0
      ;;
    d)
      fasta_path="$OPTARG"
      ;;
    d)
      if [[ $OPTARG = -* ]]; then
          >&2 echo "   ERROR: -${opt} option requires a value"
          exit 22
      fi
      extension="$OPTARG"
      ;;
    m)
      cov_mode="$OPTARG"
      #0: bidirectional, alignment covers at least coverage% of max(query, target)
      #1: alignment covers at least coverage% of target sequence
      #2: alignment covers at least coverage% of query sequence
      # 0 is OK for prokaryotes
      # for eukaryota better 2: to bring isoforms together starting from the longest isoform as query

      # from userguide:
      # With --cov-mode 0 -c [0.0,1.0] only sequences are clustered that have a sequence length overlap greater than X% of the longer of the two sequences. This coverage mode should be used to cluster full length protein sequences. Default --cluster-mode is the greedy set cover.
      # With --cov-mode 1 -c [0.0,1.0] (target-cov mode) only sequences are clustered that have a sequence length overlap greater than X% of the target sequence. The target cov mode can be used to cluster protein fragments". Default --cluster-mode is the greedy incremental clustering (by length).
      # With --cov-mode 2 -c [0.0,1.0] (query-cov mode) only sequences are clustered that have a sequence length overlap greater than X% of the query sequence. The query coverage mode can be used while searching e.g. to assure a certain level of coverage.
      # MMSeqs2 automatically picks the optimal clustering strategy based on the coverage mode (--cov-mode 0 = set cover, --cov-mode 1,2 = greedy incremental).
      ;;
    c)
      coverage="$OPTARG"
      #coverage: alignment covers at least __% of target (for cov-mode 1) or query (cov-mode 2) or both (cov-mode 0)
      ;;
    p)
      min_seq_id="$OPTARG"
      #min_seq_id: only matches above sequence identity __%
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1
      ;;
  esac
done
>&2 echo " .-- BEGUN $(date) --."

if [ -z "$fasta_path" ] || [ -z "$cov_mode" ] || [ -z "$coverage" ] || [ -z "$min_seq_id" ]; then
  >&2 echo "Options -d, -m, -c, -p are all required.\n"
  >&2 echo -e "\n$usage"
  exit 22
fi

if [ -d "$fasta_path" ]; then
  >&2 echo " |-- fasta_dir: '$fasta_path/' clustering fasta files with extension '$extension'"
else
  >&2 echo "ERROR: no such directory '$fasta_path'"
  exit 1
fi

#MAIN
setname=$(basename "$fasta_path")
combinedfasta="$BIGFASTAPATH/$setname.FASTA"

if [ -s "$combinedfasta" ]; then
  >&2 echo " |-- NOTICE: combined fasta file already exists and will be re-used:"
  >&2 echo " |-- $(ls -lh $combinedfasta)"
else
  >&2 echo " |-- Creating '$combinedfasta' concatenating files from fasta_dir, please wait..."
  sbatch --wait -J "creating combined fasta for $fasta_path" -o ${combinedfasta}.out --parsable --mem=100M --time=2:00:00 -p short --ntasks=1 --cpus-per-task=1 --nodes=1 --mail-type=FAIL,END,REQUEUE,INVALID_DEPEND,ARRAY_TASKS --mail-user=$USER --wrap="find $fasta_path -name '*${extension}' -print0 | xargs -0 cat >$combinedfasta"
  sleep 66
  >&2 echo " |-- Created:"
  >&2 echo " |-- $(ls -lh $combinedfasta)"
  #approx 40m for 120k ecoli proteome fasta files
fi
if [ -s "${combinedfasta}.idx" ]; then
  >&2 echo " |-- NOTICE: combined fasta file index already exists and will be re-used"
else
  >&2 echo " |-- Launching indexing in background:"
  sbatch -o ${combinedfasta}.idx.out -J "indexing of $(basename ${combinedfasta})" --mem=20G --time=4:00:00 -p short --ntasks=16 --cpus-per-task=1 --mail-type=FAIL,END,REQUEUE,INVALID_DEPEND,ARRAY_TASKS --mail-user=$USER --wrap="indexer.py -t 16 -v -e '>' -i '^(.+\|.*)$' -r -f ${combinedfasta} > ${combinedfasta}.idx"
fi

mmseqs_output_path="results_${fasta_path}_c${coverage}_p${min_seq_id}_cm${cov_mode}"
>&2 echo " |-- Clusters will be written under '$mmseqs_output_path'"
mkdir -pv $mmseqs_output_path

#main arguments:
command="$MMSEQS_EXEC easy-linclust $combinedfasta $mmseqs_output_path/Specie_protein $mmseqs_output_path/tmp"
#options:
command="$command --threads $threads --cov-mode $cov_mode --seq-id-mode $seq_id_mode --min-seq-id $min_seq_id -c $coverage --cluster-mode 0 --createdb-mode 1" #createdbmode: softlink fasta

>&2 echo " |-- Launching '$command'"
echo $command >$mmseqs_output_path/command
echo -e "mem\t$mem\nthreads\t$threads\ntime\t$time" >$mmseqs_output_path/slurm_params
sbatch --mem=$mem --time=$time --ntasks=1 --nodes=1 --mail-type=FAIL,END,REQUEUE,INVALID_DEPEND,ARRAY_TASKS --mail-user=$USER -o $mmseqs_output_path/slurm_${fasta_path}.out -J "easylinclust:$fasta_path" --cpus-per-task=$threads --wrap="$command; rm $mmseqs_output_path/Specie_protein_all_seqs.fasta $mmseqs_output_path/Specie_protein_rep_seq.fasta"
>&2 echo " '-- ENDED $(date) --'"
