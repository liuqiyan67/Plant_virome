#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: siRNA_analysis_preprocessing.sh
#
# PURPOSE:
# This script documents the workflow for preparing small RNA (sRNA) sequencing
# data for downstream analysis.
#
# The workflow consists of:
# 1. Downloading sRNA datasets from the SRA database.
# 2. Converting SRA files to FASTQ format.
# 3. Trimming adapters and filtering reads by length (18-30 nt) using Trimmomatic.
#
# INPUT:
#   - A list of SRA accession IDs
#   - Adapter sequence file (adapters.fa)
#
# OUTPUT:
#   - Cleaned, filtered FASTQ files ready for mapping
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Input directory containing downloaded SRA files
SRA_DIR="/path/to/sra_files"
# Output directory for raw FASTQ files
FASTQ_DIR="/path/to/raw_fastq"
# Output directory for cleaned FASTQ files
CLEAN_DIR="/path/to/clean_fastq"
# Path to Trimmomatic JAR file
TRIMMOMATIC_JAR="/path/to/trimmomatic.jar"
# Path to adapter sequences for clipping
ADAPTERS="/path/to/adapters.fa"

mkdir -p "${FASTQ_DIR}" "${CLEAN_DIR}"


# =============================================================================
# Step 1: Convert SRA to FASTQ
# =============================================================================
echo "Step 1: Converting SRA files to FASTQ format..."

# Loop through all .sra files in the input directory
for SRA_FILE in "${SRA_DIR}"/*.sra; do
    
    ACCESSION=$(basename "${SRA_FILE}" .sra)
    echo "Processing: ${ACCESSION}"
    
    # Use fasterq-dump to extract FASTQ data
    # --split-files: Write separate files for paired reads (though sRNA is usually single-end)
    fasterq-dump "${SRA_FILE}" \
                 --outdir "${FASTQ_DIR}" \
                 --threads 4
                 
done


# =============================================================================
# Step 2: Adapter Trimming and Length Filtering
# =============================================================================
echo "Step 2: Removing adapters and filtering for sRNA length (18-30 nt)..."

# Loop through all extracted FASTQ files
for FASTQ_FILE in "${FASTQ_DIR}"/*.fastq; do
    
    BASENAME=$(basename "${FASTQ_FILE}" .fastq)
    CLEAN_OUTPUT="${CLEAN_DIR}/${BASENAME}_clean.fq"
    LOG_FILE="${CLEAN_DIR}/${BASENAME}.log"
    
    echo "Trimming: ${BASENAME}"
    
    # Run Trimmomatic
    # SE: Single-end mode (typical for sRNA)
    # ILLUMINACLIP: Remove adapters (seed mismatches:2, palindrome clip:20, simple clip:7)
    # MINLEN:18 / MAXLEN:30: Keep only reads within the typical sRNA size range
    java -jar "${TRIMMOMATIC_JAR}" SE -phred33 \
         "${FASTQ_FILE}" \
         "${CLEAN_OUTPUT}" \
         ILLUMINACLIP:"${ADAPTERS}":2:20:7 \
         MINLEN:18 MAXLEN:30 \
         -threads 4 \
         2> "${LOG_FILE}"
         
    # Optional: Remove raw FASTQ to save space
    # rm "${FASTQ_FILE}"

done

echo "--- Preprocessing of sRNA data complete. Ready for mapping. ---"