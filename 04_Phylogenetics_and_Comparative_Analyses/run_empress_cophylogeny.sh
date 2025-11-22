#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_empress_cophylogeny.sh
#
# PURPOSE:
# This script documents the workflow for detailed host-virus cophylogeny analysis
# using the eMPRess tool.
#
# This analysis was specifically applied to a subset of viral clades to:
# 1. Test for significant non-random associations between virus and host phylogenies.
# 2. Visualize the topological congruence (tanglegrams).
# 3. Infer specific cophylogenetic events (Cospeciation vs. Host Switching) under
#    different event cost models.
#
# INPUT:
#   - Host tree file (*_host.nwk)
#   - Viral tree file (*.nwk)
#   - Mapping file (*.mapping) linking virus tips to host tips
#
# OUTPUT:
#   - Reconciliation CSVs (event counts)
#   - Significance test PDFs (p-values)
#   - Cost region plots
#   - Tanglegrams
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Path to the eMPRess command-line tool
EMPRESS_CLI="python empress_cli.py"
# Directory containing input tree and mapping files for cophylogeny pairs
INPUT_DIR="/path/to/cophylogeny_pairs"
# Output directory for results
OUTPUT_DIR="/path/to/empress_results"

mkdir -p "${OUTPUT_DIR}"

# Define the three cost scenarios used to test robustness
# Format: -d <dup_cost> -t <trans_cost> -l <loss_cost>
COST_LOW_T="-d 3 -t 2 -l 1"    # Scenario 1: Penalizes duplication more
COST_HIGH_T="-d 2 -t 4 -l 1"   # Scenario 2: Penalizes transfer more
COST_DEFAULT="-d 1 -t 2 -l 1"  # Scenario 3: Balanced/Default costs

# List of specific host-virus pairs selected for this analysis
PAIRS_LIST="cophylogeny_pairs_list.txt"


# =============================================================================
# Main Analysis Loop
# =============================================================================
echo "Starting eMPRess cophylogeny analysis..."

while read -r BASE_NAME; do
    
    echo "Processing pair: ${BASE_NAME}"
    
    # Define input files for this specific pair
    HOST_TREE="${INPUT_DIR}/${BASE_NAME}_host.nwk"
    VIRAL_TREE="${INPUT_DIR}/${BASE_NAME}.nwk"
    MAP_FILE="${INPUT_DIR}/${BASE_NAME}_map.mapping"
    
    
    # --- Step 1: Cost Regions Analysis ---
    # Analyze how event costs affect the reconciliation solution space
    echo "  -> Running cost regions analysis..."
    ${EMPRESS_CLI} cost-regions \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        --outfile "${OUTPUT_DIR}/${BASE_NAME}_cost_regions.pdf"
        
        
    # --- Step 2: Reconciliation (under multiple cost scenarios) ---
    # Infer the most parsimonious cophylogenetic history
    echo "  -> Running reconciliation..."
    
    # Scenario 1: Low Transfer Cost
    ${EMPRESS_CLI} reconcile \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        ${COST_LOW_T} \
        --csv "${OUTPUT_DIR}/${BASE_NAME}_recon_LowT.csv"
        
    # Scenario 2: High Transfer Cost
    ${EMPRESS_CLI} reconcile \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        ${COST_HIGH_T} \
        --csv "${OUTPUT_DIR}/${BASE_NAME}_recon_HighT.csv"
        
    # Scenario 3: Default Cost
    ${EMPRESS_CLI} reconcile \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        ${COST_DEFAULT} \
        --csv "${OUTPUT_DIR}/${BASE_NAME}_recon_Default.csv"
        
        
    # --- Step 3: Significance Testing (p-value) ---
    # Test if the host and virus trees are more congruent than expected by chance
    echo "  -> Running p-value test..."
    ${EMPRESS_CLI} p-value \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        ${COST_DEFAULT} \
        --n-samples 1000 \
        --outfile "${OUTPUT_DIR}/${BASE_NAME}_pvalue.pdf"
        
        
    # --- Step 4: Visualization (Tanglegram) ---
    # Generate a visual comparison of the two trees and their links
    echo "  -> Generating tanglegram..."
    ${EMPRESS_CLI} tanglegram \
        "${HOST_TREE}" "${VIRAL_TREE}" "${MAP_FILE}" \
        --outfile "${OUTPUT_DIR}/${BASE_NAME}_tanglegram.pdf"

done < "${PAIRS_LIST}"

echo "--- Cophylogeny analysis complete. ---"