library(ggtree)
library(tidyverse)
library(treeio)
library(ape)
library(caper)
library(furrr)


if (future::supportsMulticore()){
  
  future::plan(multicore, workers=4)
  
} else {
  
  future::plan(multisession, workers=4)
  
}

SNP_cluster_tree_dat <- 
  read_tsv('output/SNP_cluster_tree_data.tsv') %>% 
  dplyr::select(rep_genome, everything()) # need taxa column 1st for ggtree


tr <- read.raxml('RAxML_bipartitionsBranchLabels.SNP_reps')
trR <- read.tree('ROOT_DIG1.rooted.tree')

# make a ggtree object
ggtr <- ggtree(trR) 

# attach metadata to the tree
ggtr <- ggtr %<+% SNP_cluster_tree_dat

# plot tree with node labels to ID interesting clades
ggtr + geom_nodelab(aes(label=node), size=3, nudge_x = -.00002) +geom_text2(aes(subset=grepl('SX', label), label=label), nudge_x = .00008)

# node 66 defines a broad clade of interest
# node 67 defines specific clade of interest



# most recent common ancestor (node number)
MRCA <- getMRCA(phy = trR, c('SX244', 'SX245'))
# only SX244 and 245 in this clade
clade.members(x = MRCA, phy =trR, tip.labels = T)
# move down the tree by two nodes to capture some more SNP clusters
reps <- clade.members(x = MRCA - 2, phy =trR, tip.labels = T)
reps <- reps[!grepl('SX', reps)]
###

### Tree showing SNP clusters selected for further analysis
ggtr +
  geom_highlight(node=MRCA - 2) +
  geom_text2(aes(label=label,subset=grepl('SX', label)), nudge_x = .00006, size=3) +
  geom_cladelabel(node=MRCA - 2, label = 'selected', offset = .00009) + 
  geom_tippoint(aes(size=num_isolates, color=hosts ))
###

SNP_cluster_tree_dat %>% filter(asm_acc %in% reps) %>% summarise(TOT=sum(num_isolates))

PDS_accs <- SNP_cluster_tree_dat %>% filter(asm_acc %in% reps) %>% pull(PDS_acc)
# usethis::use_directory('close_relative')

close_relatives_meta <- 
  hberg_meta %>% 
  filter(PDS_acc %in% PDS_accs) %>% 
  make_dest_paths('fna','assemblies/') %>% 
  make_ftp_paths('gbk_assembly_summary.tsv') %>% 
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
  write_tsv('close_relatives_ppang.tsv', col_names = FALSE)

