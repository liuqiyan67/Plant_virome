#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_bayestraits_batch.sh
#
# PURPOSE:
# This script automates the execution of BayesTraits V4 for multiple phylogenetic
# trees. It reads a list of input files (tree + data) and runs a Reversible-Jump
# MCMC (RJ-MCMC) analysis for each pair.
#
# It demonstrates how to:
# 1. Iterate through a list of analyses.
# 2. Dynamically create a command file for each run (setting the correct LogFile).
# 3. Execute BayesTraits.
#
# INPUT:
#   - A list file (bayestraits_analysis_list.txt) where each line contains:
#     <tree_file.nexus> <data_file.txt>
#   - A template command file (bayestraits_commands_template.txt)
#
# OUTPUT:
#   - BayesTraits log files (*.log.txt) and schedule files (*.Schedule.txt)
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Path to the BayesTraits V4 executable
BAYESTRAITS_EXE="./BayesTraitsV4"
# File containing the list of analyses (Tree_File Data_File)
ANALYSIS_LIST="bayestraits_analysis_list.txt"
# Template file with common BayesTraits settings
CMD_TEMPLATE="bayestraits_commands_template.txt"
# Output directory for logs
LOG_DIR="bayestraits_logs"

mkdir -p "${LOG_DIR}"

# Check if the analysis list exists
if [ ! -f "${ANALYSIS_LIST}" ]; then
    echo "Error: Analysis list file '${ANALYSIS_LIST}' not found."
    exit 1
fi


# =============================================================================
# Step 1: Batch Processing Loop
# =============================================================================
echo "Starting batch BayesTraits analysis..."

# Read the analysis list line by line
while read -r TREE_FILE DATA_FILE; do
    
    # skip empty lines
    [[ -z "$TREE_FILE" ]] && continue
    
    # Create a base name for the output log
    BASENAME=$(basename "${TREE_FILE}" .nexus)
    LOG_FILE="${LOG_DIR}/${BASENAME}.log.txt"
    
    # Create a specific command file for this run
    # We take the template and replace the placeholder or append the LogFile name
    # Note: The template provided in the prompt already has "LogFile some_log_file...",
    # so we will use 'sed' to replace that line dynamically.
    CURRENT_CMD_FILE="${LOG_DIR}/${BASENAME}.commands.txt"
    
    # Replace the generic log file name in the template with the specific one for this run
    sed "s|LogFile some_log_file.name.log.txt|LogFile ${LOG_FILE}|" "${CMD_TEMPLATE}" > "${CURRENT_CMD_FILE}"
    
    echo "---------------------------------------------------"
    echo "Running analysis for: ${BASENAME}"
    echo "Tree: ${TREE_FILE}"
    echo "Data: ${DATA_FILE}"
    echo "Log:  ${LOG_FILE}"
    
    # Execute BayesTraits
    # The command file is passed via standard input (<)
    "${BAYESTRAITS_EXE}" "${TREE_FILE}" "${DATA_FILE}" < "${CURRENT_CMD_FILE}"
    
    # Optional: Remove the temporary command file to keep things clean
    rm "${CURRENT_CMD_FILE}"

done < "${ANALYSIS_LIST}"

echo "---------------------------------------------------"
echo "Batch analysis complete. Logs are in ${LOG_DIR}"