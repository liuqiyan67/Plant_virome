# 01_virus_discovery

this directory contains scripts documenting the initial pipeline for discovering viral contigs from raw transcriptome data and performing their initial taxonomic annotation.

## workflow overview

the overall process involves three main conceptual stages. the key commands and parameters for the most complex steps are provided in the scripts within this directory.

1.  **de novo assembly**: raw sra reads were assembled using trinity.
2.  **virus identification and filtering**: a multi-step process including:
    *   removal of host and rrna sequences.
    *   filtering of contigs by length (>1000 nt).
    *   iterative blastp-based search to identify rdrp-containing contigs.
    *   manual curation and contamination checks.
    *   final re-assembly of fragmented contigs using cap3.
3.  **taxonomic annotation**: all identified rdrp protein sequences were subjected to a formal clustering analysis using mcl to assign broad taxonomic affiliations.

## script descriptions

*   **`pipeline_virus_identification.sh`**:
    this is a documentation script that outlines the sequence of commands used for the initial virus discovery pipeline, from assembly to the final set of viral contigs. it is intended as a methodological guide and is not for direct execution.

*   **`run_mcl_for_taxonomy.sh`**:
    documents the specific workflow for performing the all-vs-all diamond blastp search and subsequent mcl clustering to assign taxonomic groups to the identified rdrp sequences.
