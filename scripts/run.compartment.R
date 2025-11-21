###################################################
# Dependencies 
###################################################
library(here)
here::i_am('scripts/run.compartment.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    # source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(future)
})

###################################################
# Load data + Set up comparisons to compute 
###################################################
individual.matrix.files <- 
    list_mcool_files(pattern='.hg38.mapq_30.1000.mcool') %>%
    filter(!isMerged) %>% 
    mutate(condition=glue('{Edit}.{Celltype}.{Genotype}')) %>% 
    mutate(replicate=glue('{CloneID}.{Celltype}.{Genotype}')) %>%
    cross_join(tibble(resolution=c(10, 25, 50, 100) * 1e3)) %>% 
    nest(data=-c(condition, resolution))

###################################################
# Generate differential compartment results
###################################################
