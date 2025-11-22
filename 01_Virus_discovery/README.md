# Directory 01: Virus Discovery and Preprocessing

This directory contains the scripts and commands documenting the initial pipeline for discovering viral contigs from raw transcriptome data and performing their initial taxonomic annotation.

## Workflow Overview

The overall process described in the main manuscript's Methods section can be broken down into the following conceptual stages. The key commands and parameters for the most complex stages are provided in the scripts within this directory.

##1. De Novo Assembly:##
Raw SRA reads were assembled using Trinity. This is a standard procedure, and an example command is provided within 'pipeline_virus_identification.sh'.

##2. Virus Identification and Filtering:##
This multi-step process is documented in the script 'pipeline_virus_identification.sh'. It includes:
#   Removal of host and rRNA sequences.
#   Filtering of contigs by length (>1000 nt).
#   An iterative BLASTp-based search to identify RdRP-containing contigs.
#   Manual curation and contamination checks.
#   Final re-assembly of fragmented contigs using CAP3.

##3. Taxonomic Annotation via MCL Clustering:##
All identified RdRp protein sequences (from this study and reference databases) were then subjected to a formal clustering analysis to assign broad taxonomic affiliations. The complete workflow for this is documented in the script '06_run_mcl_for_taxonomy.sh'.

## Script Details

#   ##'pipeline_virus_identification.sh'##: This is a ##documentation script## that outlines the sequence of commands used for the initial virus discovery pipeline, from assembly to the final set of viral contigs. It is not intended for direct execution without modification.

#   ##'run_mcl_for_taxonomy.sh'##: This script documents the specific workflow for performing the all-vs-all DIAMOND search and subsequent MCL clustering to assign taxonomic groups to the identified RdRp sequences.
