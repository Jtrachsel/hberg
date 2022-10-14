library(pdtools)
library(usethis)
library(tidyverse)
library(furrr)




if (future::supportsMulticore()){
    
    future::plan(multicore, workers=4)
    
  } else {
    
    future::plan(multisession, workers=4)
    
  }


# SX245 = JF6X01.0523
# SX244 = JF6X01.0590

# currently lacking a location for the genomes of interest

use_directory('data')
use_directory('output')


# download the most recent salmonella metadata from ncbi pathogens
dl <- download_most_recent_complete('Salmonella', folder_prefix = 'data/')


# download the genbank assembly summary file
download_gbk_assembly_summary(organism = 'Salmonella_enterica', destfile = 'gbk_assembly_summary.tsv')

# read in the metadata and join the SNP cluster designations
# this list files line is here because we could have multiple versions
# of the metadata if we are re-running these scripts over time
PDG_files <- list.files('data', 'PDG',full.names = T) %>% sort(decreasing = T)
PDD_meta <- read_tsv(PDG_files[2]) %>% 
  left_join(read_tsv(PDG_files[1]))




# Serotypes are automatically computed for salmonella isolates
# this block extracts the serotype from the'computed_types' column
PDD_meta <- 
  PDD_meta %>% 
  mutate(Serotype=sub('serotype=(.*),antigen_formula=(.*)', '\\1', computed_types), 
         antigen_formula=sub('serotype=(.*),antigen_formula=(.*)', '\\2', computed_types))


# this block filters the metadata to only Hberg genomes
# then selects the SNP clusters these genomes belong to
# and returns the unique SNP cluster accessions
# the result is every SNP cluster containing a genome with the Hberg designation
Hberg_PDS <- 
  PDD_meta %>% 
  filter(Serotype == 'Heidelberg') %>%
  filter(!is.na(PDS_acc)) %>% 
  pull(PDS_acc) %>% 
  unique()


# this filters all Salmonella to only those SNP clusters that contain 
hberg_adjacent <- PDD_meta %>% filter(PDS_acc %in% Hberg_PDS)

# check the serotypes in this collection
hberg_adjacent %>% pull(Serotype) %>% unique()

# look at different serotype designations within SNP clusters
# what proportion are Heidelberg?
WACK_PDS <- 
  hberg_adjacent %>% 
  group_by(PDS_acc, Serotype) %>%
  tally() %>% 
  summarise(non_hberg=sum(n[Serotype != 'Heidelberg']),
            hberg=sum(n[Serotype == 'Heidelberg']), 
            P_hberg=hberg/sum(hberg, non_hberg)) %>% arrange((P_hberg)) %>% 
  filter(P_hberg < 0.20) %>%  # keep all SNP clusters that are at least 20% hberg serotypes
  pull(PDS_acc)

# problem clusters?
hberg_adjacent %>%
  filter(PDS_acc %in% WACK_PDS) %>%
  pull(Serotype)


hberg_meta <- 
  hberg_adjacent %>%
  filter(!(PDS_acc %in% WACK_PDS)) %>%  # remove PDS with low % hberg
  filter(asm_acc != 'NULL')             # remove genomes without assemblies



# extract simplified metadata from existing columns
host_dat <- hberg_meta %>% extract_consensus_ag_species()
year_dat <- hberg_meta %>% extract_earliest_year()
state_dat <- hberg_meta %>% extract_state()
agency_dat <- hberg_meta %>% extract_collection_agency()
country_dat <- hberg_meta %>% extract_country()

# join these data back to the original metadata
# we should look over these dataframes to make sure they make sense
hberg_meta <- 
  hberg_meta %>%
  left_join(host_dat) %>% 
  left_join(year_dat) %>% 
  left_join(state_dat) %>% 
  left_join(agency_dat) %>% 
  left_join(country_dat)

hberg_meta %>% write_tsv('./output/01_all_hberg_metadata.tsv')
# make a summary of each SNP cluster and select a representative genome
set.seed(5)
SNP_cluster_summary <- 
  hberg_meta %>%
  group_by(PDS_acc) %>% 
  summarise(rep_genome=sample(asm_acc, 1), 
            hosts=paste(sort(unique(ag_match)), collapse = ';'), 
            years=paste(sort(unique(Year)), collapse = ';'), 
            states=paste(sort(unique(State[State !=''])), collapse = ';'),
            countries=paste(sort(unique(country[country !=''])), collapse = ';'),
            agencies=paste(sort(unique(collection_agency)), collapse = ';'), 
            num_isolates=n(),
            asm_level=asm_level[asm_acc == rep_genome],
            PFGE_patterns=paste(sort(unique(PFGE_PrimaryEnzyme_pattern)), collapse = ';'),
            .groups = 'drop')


SNP_cluster_summary$asm_level %>% unique()
#
SNP_cluster_tree_dat <- 
  SNP_cluster_summary %>%
  mutate(asm_acc=rep_genome) %>% # these functions require a column named 'asm_acc'
  make_ftp_paths('gbk_assembly_summary.tsv') %>% 
  make_download_urls('fna') %>% 
  make_dest_paths('fna', 'assemblies/') %>% 
  download_genomes(type='fna', PARALLEL = TRUE)


write_tsv(SNP_cluster_tree_dat, 'output/SNP_cluster_tree_data.tsv')
read_tsv('output/SNP_cluster_tree_data.tsv')

# decompress .gz files
gzfiles <- list.files('./assemblies/', pattern = '*gz', full.names = TRUE)
gunzip_results <- furrr::future_map(.x = gzfiles, .f = ~R.utils::gunzip(.x))
# any that dont gunzip properly?
gzfiles <- list.files('./assemblies/', pattern = '*gz')
print(gzfiles)

# identify 'complete' genomes
complete_genomes_paths <- 
  SNP_cluster_tree_dat %>% 
  filter(asm_level == 'Complete Genome') %>%
  mutate(dest=sub('.gz', '', fna_dest)) %>%
  pull(dest)

# add reference 'SL476'
complete_genomes_paths <- c(complete_genomes_paths, list.files('assemblies', pattern = 'SL476', full.names = T))

# all others are incomplete
incomplete_genomes_paths <- 
  SNP_cluster_tree_dat %>% 
  filter(asm_level != 'Complete Genome') %>%
  mutate(dest=sub('.gz', '', fna_dest)) %>%
  pull(dest)

# add in SX genomes
incomplete_genomes_paths <- c(incomplete_genomes_paths,list.files('assemblies', pattern = 'SX', full.names = T))

# build a file needed by ppanggolin to calculate a pangenome
build_ppanggolin_file_fastas(complete_genome_paths = complete_genomes_paths, 
                             incomplete_genome_paths = incomplete_genomes_paths) %>% 
  write_tsv('./output/hberg_ppanggolin_file.tsv', col_names = FALSE)
