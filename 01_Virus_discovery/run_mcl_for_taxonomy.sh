#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_mcl_for_taxonomy.sh
#
# PURPOSE:
# This script documents the workflow for taxonomic annotation of RdRp
# sequences using a DIAMOND+MCL clustering approach. It takes a FASTA file
# of RdRp proteins as input and produces a final MCL cluster file.
#
# USAGE:
# This script is for documentation. To run, replace the placeholder
# variables and execute the commands sequentially.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define variables ---
# These would be set before running the commands.

# Input FASTA file containing all RdRp sequences (newly found + reference)
INPUT_FASTA="tsa_eve_kclu50_addref.fa"

# Inflation factor for the MCL algorithm
INFLATION_FACTOR="2.0" # As per the Methods section

# Base name for output files, derived from the input file
BASENAME=$(basename "$INPUT_FASTA" .fa)

echo "--- Starting MCL clustering for: ${BASENAME} ---"
echo "--- Using Inflation Factor: ${INFLATION_FACTOR} ---"


# --- Step 1: Create DIAMOND database ---
echo "Step 1: Creating DIAMOND database..."
diamond makedb --in "${INPUT_FASTA}" -d "${BASENAME}.dmnd"


# --- Step 2: Perform all-vs-all DIAMOND BLASTp search ---
echo "Step 2: Running all-vs-all DIAMOND BLASTp..."
diamond blastp \
    --db "${BASENAME}.dmnd" \
    --query "${INPUT_FASTA}" \
    --out "${BASENAME}.blastp.out" \
    --max-target-seqs 1000000 \
    --evalue 1e-5 \
    --outfmt 6 qseqid sseqid evalue \
    --threads 150

# Check if BLASTp was successful
if [ ! -s "${BASENAME}.blastp.out" ]; then
    echo "ERROR: DIAMOND BLASTp failed or produced an empty output file. Exiting."
    exit 1
fi


# --- Step 3: Prepare input for MCL ---
echo "Step 3: Formatting BLASTp output for MCL..."
# We only need the query, subject, and e-value for mcxload
cut -f 1,2,3 "${BASENAME}.blastp.out" > "${BASENAME}.abc"


# --- Step 4: Run mcxload to create the MCL matrix ---
echo "Step 4: Running mcxload to prepare MCL matrix..."
mcxload \
    -abc "${BASENAME}.abc" \
    --stream-mirror \
    --stream-neg-log10 \
    -stream-tf 'ceil(200)' \
    -o "${BASENAME}.mci" \
    -write-tab "${BASENAME}.tab"


# --- Step 5: Run MCL for final clustering ---
echo "Step 5: Running MCL for final clustering..."
mcl "${BASENAME}.mci" \
    -use-tab "${BASENAME}.tab" \
    -I "${INFLATION_FACTOR}" \
    -o "${BASENAME}.mcl"

# Check if MCL was successful
if [ ! -s "${BASENAME}.mcl" ]; then
    echo "ERROR: MCL failed or produced an empty output file. Exiting."
    exit 1
fi

echo "--- MCL clustering workflow complete! ---"
echo "Final cluster file: ${BASENAME}.mcl"