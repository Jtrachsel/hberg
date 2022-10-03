#!/bin/bash
# must be run with "bash -i" because of 'conda activate'
set -e

conda activate ppanggolin



ppanggolin annotate --fasta ./output/hberg_ppanggolin_file.tsv --cpu 40 -o SNP_reps_pan -f
ppanggolin cluster -p SNP_reps_pan/pangenome.h5 --cpu 40
ppanggolin graph -p SNP_reps_pan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p SNP_reps_pan/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p SNP_reps_pan/pangenome.h5 -o SNP_reps_pan/MSA --cpu 40 -f
ppanggolin rgp -p SNP_reps_pan/pangenome.h5 -c 16
ppanggolin module -p SNP_reps_pan/pangenome.h5 -c 16
ppanggolin spot -p SNP_reps_pan/pangenome.h5 -c 16
ppanggolin write -o SNP_reps_pan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p SNP_reps_pan/pangenome.h5 -f
ppanggolin fasta -p ./SNP_reps_pan/pangenome.h5 --output SNP_reps_pan/REP_PROTS --prot_families all -f



raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n SNP_reps_pan -s ./SNP_reps_pan/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7


# installed root_digger to find root of tree

#rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40
rd --msa ./SNP_reps_pan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.SNP_reps_pan --prefix ROOT_DIG1 --threads 40 --exhaustive --early-stop


### TRY IQ TREE, has built in rooting options
# infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
#iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO

# infer the locus trees
#iqtree -S ALN_DIR --prefix loci -T AUTO

# compute concordance factors
#iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10