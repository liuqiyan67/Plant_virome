# 03_diversity_and_network_analyses

this directory contains scripts for quantifying viral diversity, assessing sampling saturation, and reconstructing global sequence similarity networks.

## script descriptions

*   **`run_mmseqs_clustering.sh`**:
    documents the command used to cluster viral rdrp sequences at various amino acid identity (aai) thresholds (0.9, 0.7, 0.5, 0.3) using mmseqs2.

*   **`calculate_ari.py`**:
    a python script that calculates the adjusted rand index (ari) to quantify the congruence between our operational aai clustering thresholds and formal ictv taxonomic ranks.

*   **`rarefaction_curves.r`**:
    performs the volume-weighted rarefaction analysis. it uses the `inext` package to generate diversity accumulation curves based on standardized sequencing units (ssus), correcting for uneven sequencing depth.

*   **`network_analysis.r`**:
    constructs and visualizes the rdrp sequence similarity network using the `igraph` package. it also calculates the global modularity score.

*   **`network_subsampling_analysis.r`**:
    performs a stratified subsampling experiment. it generates balanced sub-networks (equal nodes per phylum) to test if the observed high modularity is robust to sampling bias.

*   **`nmds_permanova_analysis.r`**:
    performs non-metric multidimensional scaling (nmds) and permanova tests to analyze the dissimilarity of viral community composition across different host groups.