library(tidyverse)



read_tsv('close_relatives_pan/WRITE/projection/GCA_001474545.1.tsv')


read_tsv('close_relatives_pan/WRITE/functional_modules.tsv')


read_tsv('close_relatives_pan/WRITE/modules_in_organisms.tsv')
read_tsv('close_relatives_pan/WRITE/modules_in_organisms.tsv') %>%pull(completion) %>%  hist()


read_tsv('close_relatives_pan/WRITE/modules_RGP_lists.tsv')
read_tsv('close_relatives_pan/WRITE/modules_summary.tsv')

# is this missing the genes present in these plastic regions?
read_tsv('close_relatives_pan/WRITE/plastic_regions.tsv')

read_tsv('close_relatives_pan/WRITE/')
