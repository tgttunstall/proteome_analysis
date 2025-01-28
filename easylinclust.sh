#!/usr/bin/bash

# MMseqs2 executable path
MMSEQS_EXEC="/hps/software/users/jlees/tanushree/MMseqs2/build/bin/mmseqs"
#MMSEQS_EXEC="/nfs/production/martin/uniprot/users/insana/proteomescomparisons/mmseqs-linux-avx2/bin/mmseqs"

# Input and output paths
INPUT_DIR="/nfs/research/martin/uniprot/research/spneumo_dataset/test_ds"
OUTPUT_DIR="/nfs/research/martin/uniprot/research/spneumo_dataset/test_ds/mmseqs_results"
COMBINED_FASTA="/nfs/research/martin/uniprot/research/spneumo_dataset/test_ds/all_proteins.faa"

# Create OUTPUT_DIR if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# MMseqs2 parameters
THREADS=8
COV_MODE=1
SEQ_ID_MODE=0
MIN_SEQ_ID=0.90
COVERAGE=0.80
CLUSTER_MODE=2

# Combine all FASTA files
cat ${INPUT_DIR}/*.faa > ${COMBINED_FASTA}

# Run MMseqs2
${MMSEQS_EXEC} easy-linclust ${COMBINED_FASTA} ${OUTPUT_DIR}/results ${OUTPUT_DIR}/tmp \
    --threads ${THREADS} \
    --cov-mode ${COV_MODE} \
    --seq-id-mode ${SEQ_ID_MODE} \
    --min-seq-id ${MIN_SEQ_ID} \
    -c ${COVERAGE} \
    --cluster-mode ${CLUSTER_MODE}

# Clean up
#rm ${COMBINED_FASTA}
#rm -rf ${OUTPUT_DIR}/tmp
