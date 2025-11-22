# 02_quality_control

this directory contains scripts for the systematic assessment of host veracity. the workflow integrates metagenomic quality metrics with independent functional evidence to classify viral clades into confidence tiers.

## script descriptions

*   **`calculate_contig_ratio.sh`**:
    parses the taxonomy of assembled contigs to calculate the proportion of plant-derived sequences in each library. this metric is central to our "high library quality" validation path.

*   **`calculate_virus_rpm.sh`**:
    maps reads back to viral contigs to calculate viral abundance (reads per million, rpm). this metric supports our "high viral abundance" validation path.

*   **`virus_qc_framwork.r`**:
    the core r script that integrates the metrics above. it generates the statistical plots (strip/box plots) comparing "trusted" vs. "other" datasets and produces the venn diagrams used to visualize the filtering thresholds.

*   **`sirna_analysis_preprocessing.sh`**:
    documents the initial processing of small rna (srna) datasets, including adapter trimming and size filtering.

*   **`sirna_analysis_mapping.sh`**:
    documents the mapping of cleaned srna reads to viral genomes using bowtie2.

*   **`sirna_analysis_summary.sh`**:
    aggregates the mapping results to calculate key sirna signature metrics: read length distribution (21-22 nt), genome coverage, and strand bias.