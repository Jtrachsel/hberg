#!/bin/bash
# must be run with "bash -i" because of 'conda activate'
set -e

conda activate ppanggolin



ppanggolin annotate --fasta close_relatives_ppang.tsv --cpu 40 -o close_relatives_pan -f
ppanggolin cluster -p close_relatives_pan/pangenome.h5 --cpu 40
ppanggolin graph -p close_relatives_pan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p close_relatives_pan/pangenome.h5
ppanggolin msa --source dna --partition core --phylo -p close_relatives_pan/pangenome.h5 -o close_relatives_pan/MSA --cpu 40 -f
ppanggolin rgp -p close_relatives_pan/pangenome.h5 -c 16
ppanggolin module -p close_relatives_pan/pangenome.h5 -c 16
ppanggolin spot -p close_relatives_pan/pangenome.h5 -c 16
ppanggolin write -o close_relatives_pan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p close_relatives_pan/pangenome.h5 -f
ppanggolin fasta -p ./close_relatives_pan/pangenome.h5 --output close_relatives_pan/REP_PROTS --prot_families all -f



raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n close_relatives_pan -s ./close_relatives_pan/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7


# installed root_digger to find root of tree

#rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40
rd --msa ./close_relatives_pan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.close_relatives_pan --prefix close_rel_root --threads 40 --exhaustive --early-stop


### TRY IQ TREE, has built in rooting options
# infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
#iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO

# infer the locus trees
#iqtree -S ALN_DIR --prefix loci -T AUTO

# compute concordance factors
#iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10