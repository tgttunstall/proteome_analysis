#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH -J "BATCH_ARRAY"
#SBATCH -o "array_logs/clustal_align_%A_%a.out"
#SBATCH --array=1-6%4
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_list_file> <output_directory>"
    exit 1
fi

FASTA_LIST=$1
OUTPUT_DIR=$2

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Get the FASTA file for this job
FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FASTA_LIST")

if [ -f "$FASTA_FILE" ]; then
    BASE_NAME=$(basename "$FASTA_FILE" .fa)
    BASE_NAME=${BASE_NAME%.fasta}
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}.aln"

    echo "Processing: $FASTA_FILE"
    echo "Output file: $OUTPUT_FILE"
    
    # Run Clustal Omega
    echo "clustalo -i "$FASTA_FILE" -o "$OUTPUT_FILE" -t Protein --outfmt fa --threads $SLURM_CPUS_PER_TASK --force "
    
    if [ $? -eq 0 ]; then
        echo "Successfully aligned: $FASTA_FILE"
    else
        echo "Error aligning: $FASTA_FILE"
    fi
else
    echo "File not found: $FASTA_FILE"
fi


##########
# FASTA_DIR="~/clusters_pig"
# OUTDIR="~/clusters_pig/aligned"
# THREADS=15

# # Get list of FASTA files
# FILES=(${FASTA_DIR}/*.fa ${FASTA_DIR}/*.fasta)
# TOTAL_FILES=${#FILES[@]}
# FILES_PER_JOB=$((TOTAL_FILES / 10 + 1))

# # Calculate start and end indices for this job
# START_INDEX=$((SLURM_ARRAY_TASK_ID * FILES_PER_JOB))
# END_INDEX=$((START_INDEX + FILES_PER_JOB - 1))

# # Process files for this job
# for ((i=START_INDEX; i<=END_INDEX && i<TOTAL_FILES; i++)); do
#     fasta_file="${FILES[i]}"
#     ./fasta_aligner.sh "${fasta_file}" "${OUTDIR}" --cpus-per-task=${THREADS}
# done
