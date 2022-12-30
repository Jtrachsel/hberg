#!/bin/bash
set -m 
set -e
# needs gnu parallel and bbtools 

# DEDUPE GUIDED CONTIGS #
# ncbi pathogen detection does some guided assemblies to make sure
# they dont miss AMR genes
# this sometimes results in redundant contigs, 
# or contigs that are totaly contained within other contigs
# this removes them with bbtools dedupe
conda activate ppanggolin

mkdir -p dedupe_assemblies
parallel -j 16 "dedupe.sh in={} out=dedupe_assemblies/{/.}.fna threads=1" ::: assemblies/*fna

# mv assemblies/*DD.fna dedupe_assemblies/

