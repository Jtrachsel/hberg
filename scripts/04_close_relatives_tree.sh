#!/bin/bash
# must be run with "bash -i" because of 'conda activate'
set -e

conda activate ppanggolin



ppanggolin annotate --fasta output/close_relatives_ppang.tsv --cpu 40 -o output/close_relatives_pan -f
ppanggolin cluster -p output/close_relatives_pan/pangenome.h5 --cpu 40
ppanggolin graph -p output/close_relatives_pan/pangenome.h5 -c 16
ppanggolin partition --cpu 4 -p output/close_relatives_pan/pangenome.h5
ppanggolin msa --source protein --partition core --phylo -p output/close_relatives_pan/pangenome.h5 -o output/close_relatives_pan/MSA --cpu 40 -f
ppanggolin rgp -p output/close_relatives_pan/pangenome.h5 -c 16
ppanggolin module -p output/close_relatives_pan/pangenome.h5 -c 16
ppanggolin spot -p output/close_relatives_pan/pangenome.h5 -c 16
ppanggolin write -o output/close_relatives_pan/WRITE --Rtab --csv --projection --stats --regions --spots --modules --borders --families_tsv --spot_modules -p output/close_relatives_pan/pangenome.h5 -f
ppanggolin fasta -p ./output/close_relatives_pan/pangenome.h5 --output output/close_relatives_pan/REP_PROTS --prot_families all -f



# raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n close_relatives_pan -s ./close_relatives_pan/MSA/core_genome_alignment.aln -T 35 -x 7 -N autoMRE -p 7


# installed root_digger to find root of tree

#rd --msa ./smallpan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.small_core --prefix ROOT_DIG1 --threads 40
# rd --msa ./close_relatives_pan/MSA/core_genome_alignment.aln --tree RAxML_bipartitions.close_relatives_pan --prefix close_rel_root --threads 40 --exhaustive --early-stop


### TRY IQ TREE, has built in rooting options
# infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
#iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO

# infer the locus trees
#iqtree -S ALN_DIR --prefix loci -T AUTO

# compute concordance factors
#iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10


# RUN AMRFINDER ON REP_PROTS
amrfinder -p output/close_relatives_pan/REP_PROTS/all_protein_families.faa  --organism Salmonella --plus --output output/close_relatives_pan_prots_amrfinder.tsv
# RUN PSORT ON REP_PROTS
psortb -n -r output/ -i output/close_relatives_pan/REP_PROTS/all_protein_families.faa -o terse -v
mv output/2022*_psortb_gramneg.txt output/close_rel_psortb.txt



# IQTREE
# remove empty msa files...
find output/close_relatives_pan/MSA/msa_core_protein/ -size 0 -print -delete

# need to remove non strict core msas...
mkdir bad_msas

# find max number of genomes
MAX=$(for x in output/close_relatives_pan/MSA/msa_core_protein/*aln; do cat $x |grep '>' | wc -l; done |sort -n| tail -n 1)


for x in output/close_relatives_pan/MSA/msa_core_protein/*aln
do

  NUM=$(cat $x |grep '>' | wc -l)
  NUM_UNIQ=$(seqkit rmdup -s < $x | grep '>' |wc -l)

  # if not all genomes are represented
  # OR all sequences are identical
  # remove the msa
  if (( $NUM != $MAX )) || (( $NUM_UNIQ == 1 )); then
     mv $x bad_msas/
  fi
  
  
done



# after removing proteins not found in all genomes, and invariant protiens, 
# we are left with 772 proteins for phylogenetic inference (10-14-2022)
iqtree -s output/close_relatives_pan/MSA/msa_core_protein/ -T AUTO --verbose -m MFP -msub nuclear --merge -B 1000 

rm -r bad_msas

