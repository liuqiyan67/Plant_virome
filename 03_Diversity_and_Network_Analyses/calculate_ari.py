#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SCRIPT: calculate_ari.py

PURPOSE:
This script calculates the Adjusted Rand Index (ARI) to quantify the congruence
between operational clustering thresholds (based on AAI) and formal ICTV taxonomic
ranks (Genus, Family, Order, Phylum).

It takes a reference taxonomy file and a set of clustering results (at different
AAI thresholds and coverage modes) as input, and outputs a summary matrix of ARI scores.

INPUT:
  - Taxonomy File (-t): A tab-separated file linking Seq_ID to Genus, Family, Order, Phylum.
  - Clustering Directory (-c): A directory containing subdirectories for each coverage mode,
    which in turn contain clustering results (e.g., 'labels_id_0.90.tsv').

OUTPUT:
  - A pivot table (TSV) summarizing ARI scores for each rank and threshold.

USAGE:
  python 03_calculate_ari.py -t reference_taxonomy.tsv -c mmseqs_clusters/ -o ari_results.tsv
"""

import pandas as pd
from sklearn.metrics import adjusted_rand_score
import argparse
import os
import re

def clean_id(seq_id):
    """Cleans sequence IDs by removing trailing whitespace or punctuation."""
    if not isinstance(seq_id, str): return seq_id
    return re.sub(r'[:.\s]+$', '', seq_id)

def calculate_ari_for_mode(taxonomy_df, base_cluster_dir, mode_name, thresholds):
    """Calculates ARI for a specific coverage mode across multiple AAI thresholds."""
    print(f"Analyzing Coverage Mode: {mode_name}")
    
    cluster_dir = os.path.join(base_cluster_dir, mode_name)
    if not os.path.isdir(cluster_dir):
        print(f"  Warning: Directory not found for {mode_name}. Skipping.")
        return None

    # Start with the taxonomy dataframe
    df_merged = taxonomy_df.copy()
    
    # Merge clustering results for each threshold
    for threshold in thresholds:
        # Assumes filename format: labels_id_0.90.tsv
        cluster_file = os.path.join(cluster_dir, f'labels_id_{threshold:.2f}.tsv')
        try:
            df_cluster = pd.read_csv(cluster_file, sep='\t', names=['Seq_ID', f'Cluster_ID_{threshold:.2f}'])
            df_cluster['Seq_ID'] = df_cluster['Seq_ID'].apply(clean_id)
            
            # Merge clustering data onto taxonomy data
            df_merged = pd.merge(df_merged, df_cluster, on='Seq_ID', how='left')

        except FileNotFoundError:
            print(f"  Warning: Cluster file not found: {cluster_file}")
            df_merged[f'Cluster_ID_{threshold:.2f}'] = None
    
    # Calculate ARI for each taxonomic rank
    taxonomic_ranks = ['genus', 'family', 'order', 'phylum']
    results_list = []
    
    # Normalize column names to lowercase for matching
    df_merged.columns = [col.lower().strip() for col in df_merged.columns]

    for rank in taxonomic_ranks:
        rank_label_col = rank
        if rank_label_col not in df_merged.columns: continue
        
        for threshold in thresholds:
            cluster_label_col = f'cluster_id_{threshold:.2f}'
            if cluster_label_col not in df_merged.columns: continue

            # Create a clean subset for calculation (remove NAs and "Unclassified" entries)
            df_calc = df_merged[[rank_label_col, cluster_label_col]].copy().dropna()
            valid_rows = df_calc[~df_calc[rank_label_col].astype(str).str.lower().isin(['na', 'unclassified'])]
            
            num_valid = len(valid_rows)
            ari_score = float('nan')
            
            # Only calculate ARI if we have valid data and variability
            if num_valid > 1 and len(valid_rows[rank_label_col].unique()) > 1 and len(valid_rows[cluster_label_col].unique()) > 1:
                ari_score = adjusted_rand_score(valid_rows[rank_label_col], valid_rows[cluster_label_col])
            
            results_list.append({
                'Coverage_Mode': mode_name,
                'Taxonomic_Rank': rank.capitalize(),
                'AAI_Threshold': f'{threshold:.2f}',
                'ARI_Score': ari_score,
                'N_Sequences': num_valid
            })
            
    return pd.DataFrame(results_list)

def main(taxonomy_file, base_cluster_dir, output_file):
    print("--- Starting ARI Congruence Analysis ---")

    # 1. Load and Standardize Taxonomy Data
    try:
        df_tax = pd.read_csv(taxonomy_file, sep='\t', on_bad_lines='warn')
        # Standardize the first column name to 'Seq_ID' for merging
        df_tax.rename(columns={df_tax.columns[0]: 'Seq_ID'}, inplace=True)
        df_tax['Seq_ID'] = df_tax['Seq_ID'].apply(clean_id)
        print(f"Loaded taxonomy reference: {len(df_tax)} sequences.")
    except Exception as e:
        print(f"Error parsing taxonomy file: {e}"); return

    # 2. Define Parameters
    # Standard AAI thresholds and coverage modes used in the study
    thresholds = [0.9, 0.7, 0.5, 0.3]
    coverage_modes = ['cov_short_80'] # Defaulting to the main mode used in the paper
    # Note: Add other modes ['cov_short_50', 'cov_bidir_50', 'cov_bidir_80'] if available

    all_results = []

    # 3. Iterate and Calculate
    for mode in coverage_modes:
        mode_results_df = calculate_ari_for_mode(df_tax, base_cluster_dir, mode, thresholds)
        if mode_results_df is not None:
            all_results.append(mode_results_df)

    # 4. Report Generation
    if not all_results:
        print("No results were generated. Please check input paths."); return
        
    final_df = pd.concat(all_results, ignore_index=True)
    
    # Create a readable pivot table
    pivot_table = final_df.pivot_table(
        index='Taxonomic_Rank',
        columns=['Coverage_Mode', 'AAI_Threshold'],
        values='ARI_Score'
    ).round(4)

    print("\n--- ARI Summary Matrix ---")
    print(pivot_table.to_string())
    
    pivot_table.to_csv(output_file, sep='\t')
    print(f"\nResults saved to: {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate Adjusted Rand Index (ARI) for viral clustering.")
    parser.add_argument("-t", "--taxonomy", required=True, help="Path to the taxonomy reference file (SeqID, Genus, Family...).")
    parser.add_argument("-c", "--base_cluster_dir", required=True, help="Directory containing MMseqs2 clustering results.")
    parser.add_argument("-o", "--output", required=True, help="Output path for the ARI summary file.")
    args = parser.parse_args()
    
    main(args.taxonomy, args.base_cluster_dir, args.output)