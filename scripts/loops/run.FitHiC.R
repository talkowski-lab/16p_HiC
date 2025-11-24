###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/run.FitHiC.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'loops', 'utils.loops.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

###################################################
# Set Up
###################################################
# numCores <- length(availableWorkers()); numCores
# plan(multisession, workers=numCores)
# Individual sample metadata
sample.metadata.df <- 
    load_sample_metadata()
# FitHiC2 parameters
hyper.params.df <- 
    # Default param values
    expand_grid(
        passes=c(1),
        nOfBins=c(100),
        biasLowerBound=c(0.5),
        biasUpperBound=c(2),
        mappabilityThres=c(1),
        upperbound=c(1000000, 5000000),
        normalization=c('weight', 'NONE'),
        resolution=c(5, 10, 25, 50, 100) * 1e3
    ) %>%
    # default recommended by FitHiC
    mutate(lowerbound=2 * resolution) %>%  
    # Only look at cis loops
    cross_join(
        tribble(
            ~contact_type, ~contact_type_arg,
            # 'both',        'All',
            # 'trans'        'interOnly',
            'cis',         'intraOnly' 
        )
    ) %>% 
    add_column(method='FitHiC')

###################################################
# Generate Consensus TAD calls
###################################################
# List all groups of libraries to call ConsensusTADs() from
sample.groups.df <- 
    tibble(Sample.Group=CNV.GROUPS) %>% 
    set_up_sample_groups(use_merged=FALSE)
# sample.groups.df %>% print(n=Inf)
# hyper.params.df
# Call Consensus TADs from individual groups of libraries
sample.groups.df %>%
    run_all_loops(
        hyper.params.df=hyper.params.df,
        # force_redo=TRUE,
        chromosomes=CHROMOSOMES
    )

