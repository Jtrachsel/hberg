library(tidyverse)



read_tsv('close_relatives_pan/WRITE/projection/GCA_001474545.1.tsv')


read_tsv('close_relatives_pan/WRITE/functional_modules.tsv')


read_tsv('close_relatives_pan/WRITE/modules_in_organisms.tsv')
read_tsv('close_relatives_pan/WRITE/modules_in_organisms.tsv') %>%pull(completion) %>%  hist()


read_tsv('close_relatives_pan/WRITE/modules_RGP_lists.tsv')
read_tsv('close_relatives_pan/WRITE/modules_summary.tsv')

# is this missing the genes present in these plastic regions?
read_tsv('close_relatives_pan/WRITE/plastic_regions.tsv')

# read_tsv('close_relatives_pan/WRITE/')




gpa <- read_tsv('output/close_relatives_pan/WRITE/gene_presence_absence.Rtab')




pan_mat <- gpa %>% column_to_rownames(var='Gene') %>% as.matrix()


# pan_reps <- get_pangenome_representatives(pan_mat = pan_mat, desired_coverage = .99)
library(furrr)

if (future::supportsMulticore()){
  
  future::plan(multicore, workers=12)
  
} else {
  
  future::plan(multisession, workers=12)
  
}


# pan_reps <- 
#   tibble(SEED=1:50, 
#          REPS=future_map(.x=SEED,
#                          .options = furrr_options(seed = TRUE),
#                          .f=~get_pangenome_representatives(pan_mat, SEED = .x)))

ncol(pan_mat)
nrow(pan_mat)
pangenome_reps <-
  pick_derep_sets(pan_PA = pan_mat,
                  num_sets = 100,
                  desired_coverage = .99)


novelty_scores <- calculate_novelty(pangenome_reps)

meta <- meta %>% left_join(novelty_scores %>% dplyr::select(asm_acc, log_novelty, RANK))


