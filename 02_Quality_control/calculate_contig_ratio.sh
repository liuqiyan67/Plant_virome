#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: calculate_contig_ratio.sh
#
# PURPOSE:
# This script parses taxonomy lineage information for contigs in a library
# and calculates the proportion of sequences belonging to major eukaryotic
# groups (Land plants, Algae, Fungi, Animals, etc.).
#
# INPUT:
#   - Tab-delimited file where the last column is the taxonomy lineage string.
#
# OUTPUT:
#   - A summary file (.category) with counts and percentages for each group.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to process a single file
process_file() {
    local input_file="$1"
    # Create output filename by replacing extension
    local output_file="${input_file%.merge}.category"
    
    # Initialize counters
    local total_euk=0
    local animals=0
    local fungi=0
    local landplants=0
    local allplants=0
    local rhodophyta=0
    local glaucocystophyceae=0
    
    # Read file line by line
    while IFS=$'\t' read -r -a line; do
        # Get the taxonomy lineage from the last column
        local tax_lineage="${line[-1]}"
        
        # Check if it's Eukaryota
        if [[ "$tax_lineage" == *"cellular organisms;Eukaryota;"* ]]; then
            ((total_euk++))
            
            # Classify into specific groups
            if [[ "$tax_lineage" == *"Opisthokonta;Metazoa"* ]]; then
                ((animals++))
            elif [[ "$tax_lineage" == *"Opisthokonta;Fungi"* ]]; then
                ((fungi++))
            elif [[ "$tax_lineage" == *"Embryophyta"* ]]; then
                ((landplants++))
            fi
            
            # Collect all plant-related groups to calculate Algae later
            if [[ "$tax_lineage" == *"Viridiplantae"* ]]; then ((allplants++)); fi
            if [[ "$tax_lineage" == *"Rhodophyta"* ]]; then ((rhodophyta++)); fi
            if [[ "$tax_lineage" == *"Glaucocystophyceae"* ]]; then ((glaucocystophyceae++)); fi
        fi
    done < "$input_file"
    
    # Calculate derived groups
    # Algae = (All Viridiplantae + Rhodophyta + Glaucophyta) - Land Plants
    local archaeplastida=$((allplants + rhodophyta + glaucocystophyceae))
    local algae=$((archaeplastida - landplants))
    if ((algae < 0)); then algae=0; fi
    
    local other_euk=$((total_euk - animals - fungi - archaeplastida))
    
    # Calculate percentages and write output
    {
        echo "File: $input_file"
        echo "Total Eukaryotes: $total_euk"
        
        if ((total_euk > 0)); then
            printf "Animals: %d (%.2f%%)\n" "$animals" "$(echo "scale=4; ($animals/$total_euk)*100" | bc)"
            printf "Fungi: %d (%.2f%%)\n" "$fungi" "$(echo "scale=4; ($fungi/$total_euk)*100" | bc)"
            printf "Landplants: %d (%.2f%%)\n" "$landplants" "$(echo "scale=4; ($landplants/$total_euk)*100" | bc)"
            printf "Algae: %d (%.2f%%)\n" "$algae" "$(echo "scale=4; ($algae/$total_euk)*100" | bc)"
            printf "Other Euk: %d (%.2f%%)\n" "$other_euk" "$(echo "scale=4; ($other_euk/$total_euk)*100" | bc)"
        else
            echo "No eukaryotic sequences found."
        fi
        echo "------------------------"
    } > "$output_file"
    
    echo "Processed: $input_file"
}

# Export function for parallel execution
export -f process_file

# Find all input files (*.merge) and run in parallel
# Adjust '-j 20' based on available CPU cores
find . -name "*.merge" | parallel -j 20 process_file {}

echo "Batch processing complete."