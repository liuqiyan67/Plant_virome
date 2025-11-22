#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_dali_structural_comparison.sh
#
# PURPOSE:
# This script documents the workflow for all-vs-all protein structure comparison
# using DaliLite. This analysis is used to infer distant evolutionary relationships
# between viral Capsid Proteins (CP) or Glycoproteins (GP) that may not be
# detectable by sequence similarity alone.
#
# The workflow consists of:
# 1. Importing PDB structure files into the Dali database format.
# 2. Running an all-vs-all comparison matrix.
#
# INPUT:
#   - A directory containing PDB files (predicted by AlphaFold2 or from PDB).
#   - A list of query structures to compare.
#
# OUTPUT:
#   - A Dali DAT directory containing processed structure data.
#   - An all-vs-all distance matrix.
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Path to DaliLite executable directory
DALI_BIN="/path/to/DaliLite.v5/bin"
# Directory containing input PDB files
PDB_DIR="/path/to/pdb_structures"
# Output directory for Dali internal data files
DAT_DIR="dali_dat_files"
# List of PDB IDs to query (formatted for Dali)
QUERY_LIST="query_list.txt"
# Output directory for results
RESULTS_DIR="dali_results"

mkdir -p "${DAT_DIR}" "${RESULTS_DIR}"


# =============================================================================
# Step 1: Import PDB Structures into Dali Database
# =============================================================================
echo "Step 1: Importing PDB files..."

# Note: In a real pipeline, this would be a loop.
# Here, we show the command structure used for each file.

# Example command structure:
# ${DALI_BIN}/import.pl --pdbfile <pdb_file> --pdbid <unique_id> --dat <dat_dir>

# Examples from the study:
${DALI_BIN}/import.pl --pdbfile "${PDB_DIR}/fold_avm87615_1_structural_polyprotein_model_0.pdb" --pdbid t001 --dat "${DAT_DIR}"
${DALI_BIN}/import.pl --pdbfile "${PDB_DIR}/fold_avm87615_1_structural_polyprotein_model_1.pdb" --pdbid t002 --dat "${DAT_DIR}"
${DALI_BIN}/import.pl --pdbfile "${PDB_DIR}/fold_dad54814_1_tpa_inf_structural_protein_model_0.pdb" --pdbid t003 --dat "${DAT_DIR}"
# ... (repeated for all structures in the dataset) ...

echo "Import complete. Structures are stored in ${DAT_DIR}."


# =============================================================================
# Step 2: Run All-vs-All Structural Comparison
# =============================================================================
echo "Step 2: Running Dali all-vs-all matrix comparison..."

# Change to results directory to keep output clean
cd "${RESULTS_DIR}" || exit

# Run Dali
# --matrix: Calculate all-vs-all comparison
# --query: List of PDB IDs to compare
# --dat1: Location of the imported data files
# --clean: Clean up temporary files
nohup "${DALI_BIN}/dali.pl" \
    --matrix \
    --query "../${QUERY_LIST}" \
    --dat1 "../${DAT_DIR}/" \
    --clean \
    2> /dev/null &

echo "Dali job submitted. Check ${RESULTS_DIR} for outputs."