library(treeio)
library(ggtree)
library(tidyverse)

crtr <- read.tree('close_rel_root.rooted.tree')
read_tsv('close_relatives_ppang.tsv', col_names = F) %>% filter(grepl('SX', X1))

ggtr <- ggtree(crtr, layout = 'circular') #+ geom_tiplab() + xlim(0,.0002)


ggtr + geom_tiplab(aes(subset=grepl('SX', label)))
grep('SX', ggtr$data$label)

meta <- read_tsv('output/close_relatives_meta.tsv')

PDG <- sub('.cluster_list.tsv','',list.files('data') %>% sort(decreasing = T) %>% pluck(1))

# download NCBI snp tree for this cluster
SNP_tree_dl <- 
  tibble(
    SNP_tree_url=make_SNPtree_urls(data = meta, organism = 'Salmonella',
    PDG = PDG)
  )%>%
  make_SNP_tree_dest(data_dir = 'data') %>%
  download_SNP_trees()

meta %>%
  dplyr::select(asm_acc,PDS_acc, ag_match, Year, isolation_source, State, country, collection_agency)

meta %>%
  group_by(ag_match) %>%
  tally()

LOOK <- meta %>%
  filter(ag_match == 'Other') %>% 
  dplyr::select(outbreak,isolation_source,host, ontological_term, ag_match, Year, country, State, collection_agency)

meta %>%
  count(Year,collection_agency) %>%
  ggplot(aes(x=Year,fill=collection_agency, y=n)) +
  geom_col()


meta %>%
  count(ag_match) %>%
  ggplot(aes(x=ag_match, y=n)) +
  geom_col()

meta %>%
  group_by(Year) %>%
  tally() %>%
  ggplot(aes(x=Year, y=n)) + 
  geom_line()+
  geom_point()


meta %>% group_by(Year, ag_match) %>% tally() %>% 
  ggplot(aes(x=Year, y=n, fill=ag_match)) + geom_col()



library(usmap)
library(ggplot2)
STATE_DAT <- meta %>%
  group_by(State) %>%
  tally() %>%
  filter(!is.na(State)) %>% 
  transmute(state=State, 
            number_of_isolates=n)


plot_usmap(data = STATE_DAT, values = "number_of_isolates") + 
  scale_fill_continuous(name = "Number of isolates", label = scales::comma) + 
  scale_fill_viridis_c()+
  theme(legend.position = "right")




