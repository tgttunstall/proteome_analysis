#!/bin/bash

#SBATCH --job-name=Cpig_c0.5_p0.65_cov0
#SBATCH --array=1-12904%1000
#SBATCH --time=00:08:00
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH -o "console_clustal/Cpig_c0.5_p0.65_cov0_l%A_%a.out"   # job output file


#fasta_list="/nfs/research/martin/uniprot/research/proteome_analysis/list-thing"
#fasta_list="/nfs/research/martin/uniprot/research/proteome_analysis/clustes_pig_fasta_list.txt"
fasta_list="/nfs/research/martin/uniprot/research/proteome_analysis/Cpig_c0.5_p0.65_cov0_fasta_list.txt"

# Get the total number of lines in the FASTA list file
total_lines=$(wc -l < "$fasta_list")

sequence_type="Protein"
output_file_format="fa"
#output_dir="/nfs/research/martin/uniprot/research/proteome_analysis/test-output"
#output_dir="/nfs/research/martin/uniprot/research/proteome_analysis/results_clusters_pig
output_dir="/nfs/research/martin/uniprot/research/proteome_analysis/results_Cpig_c0.5_p0.65"

n_threads=${3:-8}
start_time_whole=$(date +%s)
date
if [ ! -d "${output_dir}" ]; then
  >&2 echo "Creating output directory: ${output_dir}"
  mkdir -p "${output_dir}"
fi

# Check if the current task ID is within the range of available files
if [ "$SLURM_ARRAY_TASK_ID" -le "$total_lines" ]; then
    # Get the FASTA file for this job
    fasta_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$fasta_list")
    # derive base_name
    base_name=$(basename ${fasta_file} .fa)
    output_file="${output_dir}/${base_name}.aln"

    
    echo "Task ID: $SLURM_ARRAY_TASK_ID is processing: $fasta_file"
    clustalo -i $fasta_file -o "${output_file}" -t "${sequence_type}" --outfmt "${output_file_format}" --threads ${n_threads} --force
    end_time_whole=$(date +%s)
    execution_time_whole=$((end_time_whole - start_time_whole))
    hours_whole=$((execution_time_whole / 3600))
    minutes_whole=$(((execution_time_whole % 3600) / 60))
    seconds_whole=$((execution_time_whole % 60))

    printf "\nAll alignments completed."
    printf "\n${fasta_count} fasta file alignments completed in %02d:%02d:%02d (HH:MM:SS)\n" ${hours_whole} ${minutes_whole} ${seconds_whole}
    date
    
    # Check if the file exists
    if [ -f "$fasta_file" ]; then
        echo "File exists and is readable."
    else
        echo "Warning: File does not exist or is not readable."
    fi
else
    echo "Task ID: $SLURM_ARRAY_TASK_ID exceeds the number of available files ($total_lines)."
fi

#SBATCH --output=test_clustal_%A_%a.out
#SBATCH --error=test_clustal_%A_%a.err
