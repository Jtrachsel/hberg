library(treeio)
library(ggtree)
library(tidyverse)

crtr <- read.tree('close_rel_root.rooted.tree')
read_tsv('close_relatives_ppang.tsv', col_names = F) %>% filter(grepl('SX', X1))

ggtr <- ggtree(crtr, layout = 'circular') #+ geom_tiplab() + xlim(0,.0002)


ggtr + geom_tiplab(aes(subset=grepl('SX', label)))
grep('SX', ggtr$data$label)

meta <- read_tsv('output/close_relatives_meta.tsv')
meta %>% dplyr::select(asm_acc,PDS_acc, ag_match, Year, isolation_source, State, country)
