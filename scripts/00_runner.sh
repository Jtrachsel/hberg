#!/bin/bash
set -e
# maybe eventually try and migrate to snakemake?

Rscript ./scripts/01_downloads.R
# outputs: 
#   directories: 
#     ./data/
#     ./output/
#     ./assemblies/
#   files:
#     - Two NCBI pathogens metadata files
#     - Genbank assembly summary file
#     - ./output/SNP_cluster_tree_data.tsv   --- representative genome + SNP cluster summary
#     - variable number of gunzipped files in ./assemblies/
#     - ./output/hberg_ppanggolin_file.tsv  --- for running SNP cluster rep pangenome RENAME
#     - ./output/01_all_hberg_metadata.tsv

bash -i ./scripts/02_SNP_cluster_tree.sh
# outputs:
#   directories:
#     ./SNP_reps_pan/

Rscript ./scripts/03_extract_close_relatives.R

bash -i ./scripts/04_close_relatives_tree.sh

Rscript ./scripts/05_plot_tree.R


