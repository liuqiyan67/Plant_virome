#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: run_phylogenetic_analysis.sh
#
# PURPOSE:
# This script documents the workflow for phylogenetic inference of viral RdRP
# sequences. It covers multiple sequence alignment, alignment trimming/refinement,
# and Maximum Likelihood (ML) tree construction.
#
# The workflow consists of:
# 1. High-accuracy local alignment using MAFFT (L-INS-i).
# 2. Automated alignment trimming using TrimAl to remove poorly aligned regions.
# 3. ML tree inference using IQ-TREE 2 with ModelFinder and SH-aLRT support.
#
# INPUT:
#   - A directory containing FASTA files of unaligned RdRP protein sequences.
#
# OUTPUT:
#   - Aligned FASTA files (*.aln)
#   - Trimmed alignment files (*.trimmed.aln)
#   - Phylogenetic tree files (*.treefile) and IQ-TREE logs.
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Input directory containing unaligned FASTA files
INPUT_DIR="/path/to/unaligned_sequences"
# Output directory for alignments
ALN_DIR="/path/to/alignments"
# Output directory for trees
TREE_DIR="/path/to/trees"

# Tools (assumed to be in PATH)
MAFFT="mafft"
TRIMAL="trimal"
IQTREE="iqtree2"

mkdir -p "${ALN_DIR}" "${TREE_DIR}"


# =============================================================================
# Step 1: Multiple Sequence Alignment
# =============================================================================
echo "Step 1: Aligning sequences using MAFFT..."

# Loop through all FASTA files in the input directory
for FASTA in "${INPUT_DIR}"/*.fa; do
    
    BASENAME=$(basename "${FASTA}" .fa)
    OUTPUT_ALN="${ALN_DIR}/${BASENAME}.aln"
    
    echo "Aligning: ${BASENAME}"
    
    # Run MAFFT
    # --localpair --maxiterate 1000: This corresponds to the L-INS-i algorithm,
    # which is highly accurate for datasets with divergent sequences.
    # --reorder: Reorders output sequences by similarity.
    # --thread -1: Auto-detect available cores.
    "${MAFFT}" --localpair --maxiterate 1000 --reorder --thread -1 \
               "${FASTA}" > "${OUTPUT_ALN}"
               
done


# =============================================================================
# Step 2: Alignment Trimming and Refinement
# =============================================================================
echo "Step 2: Trimming alignments manually..."

# =============================================================================
# Step 3: Phylogenetic Tree Inference
# =============================================================================
echo "Step 3: Inferring ML trees using IQ-TREE"

# Loop through the trimmed alignment files
for TRIMMED_ALN in "${ALN_DIR}"/*.trimmed.aln; do
    
    BASENAME=$(basename "${TRIMMED_ALN}" .trimmed.aln)
    PREFIX="${TREE_DIR}/${BASENAME}"
    
    echo "Building tree for: ${BASENAME}"
    
    # -s: Input alignment
    # -m MFP: Run ModelFinder to automatically select the best-fit substitution model
    # --alrt 1000: Perform 1000 replicates of the SH-like approximate Likelihood Ratio Test
    # -T AUTO: Auto-detect best number of CPU threads
    # --prefix: Output file prefix
    
    # Standard Run (ModelFinder + Tree Search + SH-aLRT):
    "${IQTREE}" -s "${TRIMMED_ALN}" \
                -m MFP \
                --alrt 1000 \
                -T AUTO \
                --prefix "${PREFIX}"

    # Note: For specific large datasets (e.g., Stelpaviricetes), site-heterogeneous
    # models (like -m LG+C20+F+G) were used instead of MFP.
    
done

echo "--- Phylogenetic analysis workflow complete. ---"