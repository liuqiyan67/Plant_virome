#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: siRNA_analysis_summary.sh
#
# PURPOSE:
# This script aggregates the various siRNA metrics calculated in the previous
# steps into a single, comprehensive summary report for each viral contig.
#
# It collects:
# 1. Mapped read counts
# 2. Genome coverage statistics
# 3. Strand bias (Forward/Reverse ratio)
# 4. 5' terminal nucleotide composition
# 5. Read length distribution (specifically the % of 21-22 nt reads)
#
# INPUT:
#   - Directory containing intermediate analysis files (coverage, strand counts, etc.)
#
# OUTPUT:
#   - A final summary table: siRNA_analysis_summary.txt
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Directory containing analysis results from Step 2
ANALYSIS_DIR="/path/to/analysis_results"
# Final output report file
REPORT_FILE="siRNA_analysis_summary.txt"

# Initialize the final report with a header
echo -e "Sample_ID\tVirus_Contig\tTotal_Reads\tCoverage_%\tStrand_Ratio\t21_22nt_%\t5_Prime_Bias" > "${REPORT_FILE}"

echo "Generating summary report..."


# =============================================================================
# Step 1: Aggregate Metrics for Each Sample
# =============================================================================
# Loop through each sample ID found in the analysis directory
for ID in $(ls "${ANALYSIS_DIR}"/*.idxstats.txt | xargs -n 1 basename | sed 's/.idxstats.txt//'); do
    
    echo "Summarizing: ${ID}"
    
    # 1. Get Total Mapped Reads
    TOTAL_READS=$(awk '{sum+=$3} END {print sum}' "${ANALYSIS_DIR}/${ID}.idxstats.txt")
    
    # 2. Get Genome Coverage (%)
    # Assuming a pre-calculated coverage file from the previous step
    COVERAGE=$(grep "Coverage_Rate" "${ANALYSIS_DIR}/${ID}.coverage.txt" | cut -f2)
    
    # 3. Calculate Strand Ratio (Forward / Reverse)
    FWD=$(awk '{sum+=$2} END {print sum}' "${ANALYSIS_DIR}/${ID}.forward_reads.txt")
    REV=$(awk '{sum+=$2} END {print sum}' "${ANALYSIS_DIR}/${ID}.reverse_reads.txt")
    
    if [[ "${REV}" -gt 0 ]]; then
        STRAND_RATIO=$(echo "scale=2; ${FWD} / ${REV}" | bc)
    else
        STRAND_RATIO="NA"
    fi
    
    # 4. Calculate 21-22 nt Read Percentage
    # Sum counts for lengths 21 and 22, divide by total reads
    COUNT_21=$(grep -w "21" "${ANALYSIS_DIR}/${ID}.read_length_dist.txt" | awk '{print $2}')
    COUNT_22=$(grep -w "22" "${ANALYSIS_DIR}/${ID}.read_length_dist.txt" | awk '{print $2}')
    
    # Handle missing counts (if 0)
    COUNT_21=${COUNT_21:-0}
    COUNT_22=${COUNT_22:-0}
    
    if [[ "${TOTAL_READS}" -gt 0 ]]; then
        PERCENT_21_22=$(echo "scale=2; ((${COUNT_21} + ${COUNT_22}) / ${TOTAL_READS}) * 100" | bc)
    else
        PERCENT_21_22="0"
    fi
    
    # 5. Get Major 5' Nucleotide
    # Find the base (A, C, G, T) with the highest count
    MAJOR_5P=$(sort -k2 -nr "${ANALYSIS_DIR}/${ID}.five_prime_base.txt" | head -n 1 | awk '{print $1}')
    
    
    # --- Write Row to Report ---
    echo -e "${ID}\t${ID}\t${TOTAL_READS}\t${COVERAGE}\t${STRAND_RATIO}\t${PERCENT_21_22}\t${MAJOR_5P}" >> "${REPORT_FILE}"

done

echo "--- Summary report generated: ${REPORT_FILE} ---"