###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/run.multiHiCCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'utils.multiHiCCompare.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

###################################################
# Set up all comparisons
###################################################
# Individual sample metadata
sample.metadata.df <- 
    load_sample_metadata(filter=TRUE)
# All combinations of multiHiCCompare hyper-params to test
hyper.params.df <- 
    expand_grid(
        zero.p=c(0.8),
        A.min=c(5)
    )
# List all separate sample sets + parameters to run multiHiCComapre for
comparisons.df <- 
    tribble(
        ~Sample.Group.Left, ~Sample.Group.Right,
        # ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        '16p.iN.DUP',       '16p.iN.WT',  
        '16p.iN.DEL',       '16p.iN.WT',  
        # '16p.iN.DUP',       '16p.iN.DEL', 
        '16p.NSC.DUP',      '16p.NSC.WT',
        '16p.NSC.DEL',      '16p.NSC.WT',
        # '16p.NSC.DUP',      '16p.NSC.DEL',
        # '16p.iN.DUP',       '16p.NSC.DUP',
        # '16p.iN.DEL',       '16p.NSC.DEL',
        '16p.NSC.WT',       '16p.iN.WT'
    ) %>% 
    mutate(
        across(
            starts_with('Sample.Group.'),
            ~ str_replace_all(.x, 'All', '.*'),
            .names='{.col}.Pattern'
        )
    ) %>% 
    set_up_sample_comparisons() %>%
    select(-c(ends_with('.Pattern')))

###################################################
# Generate DAC results for each comparison
###################################################
# GRanges object with Centro/Telomere regions to filter
data('hg38_cyto') 
# parallelizing params
library(BiocParallel)
# numCores <- length(availableWorkers()); numCores
numCores <- 2
# register(MulticoreParam(workers=numCores), default=TRUE)
register(MulticoreParam(workers=numCores * 2 / 4), default=TRUE)
plan(multisession,      workers=numCores * 2 / 4)
# Run multiHiCComapre on all comparisons
# 2 group comparison + no covariates -> use exact test
comparisons.df %>% 
    run_all_multiHiCCompare(    
        hyper.params.df=hyper.params.df,
        remove.regions=hg38_cyto,
        # covariates.df=NULL,
        # chromosomes=CHROMOSOMES,
        # force_redo=TRUE,
        sample_group_priority_fnc=sample_group_priority_fnc_16p
    )

