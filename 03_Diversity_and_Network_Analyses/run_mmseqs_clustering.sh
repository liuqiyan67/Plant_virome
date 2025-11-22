#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_mmseqs_clustering.sh
#
# PURPOSE:
# This script performs hierarchical clustering on a FASTA file of protein sequences
# using MMseqs2. It iterates through standard AAI thresholds (0.9, 0.7, 0.5, 0.3)
# and generates tabular output files linking each sequence to its cluster.
#
# INPUT:
#   - A FASTA file of protein sequences.
#
# OUTPUT:
#   - A directory containing TSV files for each identity threshold.
#     Each file has two columns: Sequence_ID, Cluster_Representative_ID.
#
# USAGE:
#   bash 02_run_mmseqs_clustering.sh <input_fasta> <output_directory>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: setup ---
INPUT_FASTA="$1"
OUTPUT_DIR="$2"

if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Missing arguments."
    echo "Usage: bash 02_run_mmseqs_clustering.sh <input_fasta> <output_directory>"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p "${TMP_DIR}"

# MMseqs2 Parameters (Standardized for this study)
COVERAGE_MODE=0         # Coverage of query and target
COVERAGE_THRESHOLD=0.8  # 80% coverage required
E_VALUE="1e-5"
THREADS=8

# Identity thresholds to process
THRESHOLDS=(0.90 0.70 0.50 0.30)


# --- Step 1: Create MMseqs2 Database ---
echo "Step 1: Creating sequence database..."
SEQ_DB="${TMP_DIR}/seq_db"
mmseqs createdb "${INPUT_FASTA}" "${SEQ_DB}" -v 0


# --- Step 2: Run Clustering Loop ---
echo "Step 2: Running clustering at multiple thresholds..."

for id in "${THRESHOLDS[@]}"; do
    echo "  -> Processing Identity: ${id}"
    
    CLU_DB="${TMP_DIR}/cluster_${id}"
    TSV_FILE="${OUTPUT_DIR}/labels_id_${id}.tsv"

    # Run clustering
    # --cluster-mode 1: CD-HIT-like greedy clustering
    mmseqs cluster "${SEQ_DB}" "${CLU_DB}" "${TMP_DIR}/tmp_${id}" \
        --cluster-mode 1 \
        --min-seq-id "${id}" \
        -e "${E_VALUE}" \
        --cov-mode "${COVERAGE_MODE}" \
        -c "${COVERAGE_THRESHOLD}" \
        --threads "${THREADS}" \
        -v 0

    # Convert cluster results to TSV format
    # We create a flat file linking Member -> Representative
    mmseqs createtsv "${SEQ_DB}" "${SEQ_DB}" "${CLU_DB}" "${TSV_FILE}.raw" --threads "${THREADS}" -v 0
    
    # Clean up format: Column 2 (Member) -> Column 1 (Representative)
    # Final Output Format: Sequence_ID [TAB] Cluster_ID
    awk -F'\t' '{print $2"\t"$1}' "${TSV_FILE}.raw" > "${TSV_FILE}"
    rm "${TSV_FILE}.raw"
    
    echo "     Created: ${TSV_FILE}"
done

# --- Step 3: Cleanup ---
# Uncomment the line below to remove large temporary database files after run
# rm -rf "${TMP_DIR}"

echo "--- Clustering workflow complete. Results are in ${OUTPUT_DIR} ---"