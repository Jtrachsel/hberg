library(ggtree)
library(tidyverse)
library(treeio)
library(ape)
library(caper)
library(furrr)
library(pdtools)


if (future::supportsMulticore()){
  
  future::plan(multicore, workers=4)
  
} else {
  
  future::plan(multisession, workers=4)
  
}

SNP_cluster_tree_dat <- 
  read_tsv('output/SNP_cluster_tree_data.tsv') %>% 
  dplyr::select(rep_genome, everything()) # need taxa column 1st for ggtree


trR <- read.tree('output/ROOT_DIG1.rooted.tree')



# make a ggtree object (full tree)
ggtr <- ggtree(trR) 

# attach metadata to the tree
ggtr <- ggtr %<+% SNP_cluster_tree_dat

# get most recent common ancestor (node number) for genomes of interest
MRCA <- getMRCA(phy = trR, c('SX244', 'SX245'))

# extract the clade one node up the tree from MRCA
tr_clade <- 
  extract.clade(trR, MRCA - 1 ) 

# node for the clade of interest in the subtree
MRCA_clade <- getMRCA(tr_clade, c('SX244', 'SX245'))

# make a ggtree object (clade of interest)and attach metadata to it
ggtr_clade <- 
  ggtree(tr_clade) %<+% SNP_cluster_tree_dat

# look at the node numbers on the tree
ggtr_clade +geom_nodelab(aes(label=node))

# plot full tree 
ggtr +
  # geom_nodelab(aes(label=node), size=3, nudge_x = -.00002) +
  geom_text2(aes(subset=grepl('S', label), label=label),size=3, nudge_x = .0003)+
  geom_highlight(node=MRCA - 1, extend=1 ,to.bottom = T)+
  expand_limits(x=.009)
  


# One other SNP cluster associated with SX244 and SX245
close_clusters <- clade.members(x = MRCA, phy =trR, tip.labels = T)
close_clusters <- close_clusters[!grepl('SX', close_clusters)]

### THIS ONE
p_clade <- 
  ggtr_clade + 
  # ggtree::rotate(68)+  # THIS CANT BE HARDCODED.... HOW TO CALC!?!
  # geom_nodelab(aes(label=node), size=3) +
  geom_highlight(node=MRCA_clade, extend=.000005 )+
  geom_point2(aes(subset=grepl('SX', label)),size=3, color='purple')+
  geom_tippoint(size=.3)+
  geom_tippoint(aes(size=num_isolates), alpha=.75)+
  geom_text2(aes(subset=grepl('S', label), label=label),size=3, nudge_x = .000002)+
  geom_label2(aes(subset=grepl(close_clusters, label), label=PDS_acc),size=4, nudge_x = .000011)+
  expand_limits(y=-1)+
  ggtitle('SNP cluster most closely associated with these isolates')

ggsave(p_clade, filename = 'output/SNP_cluster_tree.jpeg', width = 9, height = 7, bg='white')



# close clusters is more reasonable with 339 genomes
SNP_cluster_tree_dat %>% filter(asm_acc %in% close_clusters) %>% summarise(TOT=sum(num_isolates))

# SNP clusters to get genomes from
PDS_accs <-
  SNP_cluster_tree_dat %>%
  filter(asm_acc %in% close_clusters) %>%
  pull(PDS_acc)
# usethis::use_directory('close_relative')

# filter all hberg metadata to just those in this snp cluster
# download these genomes
close_relatives_meta <- 
  read_tsv('output/01_all_hberg_metadata.tsv') %>% 
  filter(PDS_acc %in% PDS_accs) %>% 
  make_dest_paths('fna','assemblies/') %>% 
  make_ftp_paths('data/gbk_assembly_summary.tsv') %>% 
  make_download_urls('fna') %>% 
  download_genomes(type = 'fna', PARALLEL = TRUE) %>% 
  write_tsv('output/close_relatives_meta.tsv')

SX_paths <- list.files('assemblies', 'SX', full.names = T)
# SX_dest <- paste0('close_relatives/', list.files('assemblies', 'SX'))
# file.copy(SX_paths, SX_dest)


gzfiles <- list.files('./assemblies', pattern = '*gz', full.names = TRUE)
gunzip_results <- furrr::future_map(.x = gzfiles, .f = ~R.utils::gunzip(.x))
# remove those that dont gunzip properly
gzfiles <- list.files('./assemblies', pattern = '*gz')
print(gzfiles)

# this should output the paths of the assemblies listed in the close_relatives metadata
close_rel_pattern <- paste(close_relatives_meta$asm_acc, collapse = '|')
close_relatives_paths <- 
  grep(close_rel_pattern, 
       list.files('assemblies', full.names = T),
       value = TRUE)

close_relatives_paths <- c(close_relatives_paths, SX_paths)

build_ppanggolin_file_fastas(incomplete_genome_paths =close_relatives_paths ) %>% 
  write_tsv('output/close_relatives_ppang.tsv', col_names = FALSE)

