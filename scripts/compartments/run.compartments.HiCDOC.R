###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/compartments/run.compartments.HiCDOC.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

###################################################
# Set up all hyper-params 
###################################################
options(scipen=999)
# dcHiC hyper-params
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# HiCDOC hyper-params
hyper.params.df <- 
    expand_grid(
        chrThreshold=c(100),   # Keeping chromosomes with at least 100 positions.
        repThreshold=c(0.005), # Keeping replicates filled with at least 5% non-zero interactions 
        posThreshold=c(0.01),  # Keeping positions with interactions average greater or equal to 1
        loessSampleSize=c(20000)
    )

###################################################
# Generate differential compartment results 
###################################################
# List all pairs of matrices to compare
comparisons.df <- 
    tribble(
        # ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        ~Sample.Group.Right, ~Sample.Group.Left,
        # '16p.iN.DUP',       '16p.iN.DEL', 
        # '16p.NSC.DUP',      '16p.NSC.DEL',
        # '16p.NSC.DUP',      '16p.iN.DUP',
        # '16p.NSC.DEL',      '16p.iN.DEL',
        # '16p.NSC.WT',       '16p.iN.WT',
        # '16p.iN.DUP',       '16p.iN.WT',  
        # '16p.iN.DEL',       '16p.iN.WT',  
        '16p.NSC.DUP',      '16p.NSC.WT',
        '16p.NSC.DEL',      '16p.NSC.WT'
    ) %>% 
    set_up_sample_comparisons(
        resolutions=parsed.args$resolutions,
        merging='individual'
    ) %>% 
    # This is necessary for HiCDOC input
    rowwise() %>%
    mutate(
        samples.df=
            samples.df %>%
            dplyr::rename('conditions'=Sample.Group) %>%
            group_by(conditions) %>% 
            mutate(replicates=row_number()) %>%
            ungroup() %>% 
            select(-c(resolution)) %>% 
            list()
    ) %>% 
    ungroup()
# parallelizing params
# used by calls to future_pmap() in functions below
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
# generate compartments
comparisons.df %>% 
    run_all_HiCDOC(    
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo,
        sample_group_priority_fnc=sample_group_priority_fnc_16p
    )

