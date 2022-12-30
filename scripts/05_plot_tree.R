library(treeio)
library(ggtree)
library(tidyverse)
library(pdtools)
# crtr <- read.tree('close_rel_root.rooted.tree')
crtr <- read.tree('output/close_relatives_pan/MSA/msa_core_protein.treefile')
meta <-
  read_tsv('output/close_relatives_meta.tsv') %>% 
  dplyr::select(asm_acc, everything())

# meta <- meta %>% left_join(novelty_scores %>% dplyr::select(asm_acc, log_novelty, RANK))


# ggtr <- ggtree(crtr, layout = 'circular') %<+% meta
ggtr <- ggtree(crtr) %<+% meta


hist(crtr$edge)

# ggtr$data$sra_center
ggtr +
  # geom_tippoint(aes(fill=ag_match),size=2, shape =21)+
  geom_tiplab(aes(subset=grepl('SX', label)), align = T) + 
  scale_fill_brewer(palette = 'Dark2', direction = -1) + 
  geom_point2(aes(subset=sra_center == 'USDA-NVSL-DBL'), shape=21,fill='red')+
  # geom_nodepoint(aes(subset=bootstrap > 74, color=bootstrap)) + 
  expand_limits(x=.0004) + ggtitle('')

# ggtr +
#   # geom_tippoint(aes(fill=ag_match),size=2, shape =21)+
#   geom_tiplab(aes(subset=grepl('SX', label)), align = T) + 
#   scale_fill_brewer(palette = 'Dark2', direction = -1) + 
#   geom_point2(aes(subset=!is.na(log_novelty) ),size=3, shape=21,fill='blue')+
#   geom_point2(aes(subset=sra_center == 'USDA-NVSL-DBL'), shape=21,fill='red')+
#   # geom_nodepoint(aes(subset=bootstrap > 74, color=bootstrap)) + 
#   expand_limits(x=.0004) + ggtitle('')


ggtr$data <- 
  ggtr$data %>%
  mutate(bootstrap=ifelse(as.numeric(label) > 95, as.numeric(label), NA))

ggtr$data$label




meta <- read_tsv('output/close_relatives_meta.tsv')

PDG <- sub('.cluster_list.tsv','',list.files('data', 'PDG') %>% sort(decreasing = T) %>% pluck(1))

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




