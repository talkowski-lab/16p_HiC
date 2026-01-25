###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/DifferentialContacts/run.multiHiCCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts', 'locations.R'))
    source(file.path(BASE_DIR,   'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts', 'utils.multiHiCCompare.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(optparse)
    library(BiocParallel)
    library(hictkR)
})

###################################################
# Set up all comparisons
###################################################
options(scipen=999)
# All combinations of multiHiCCompare hyper-params to test
hyper.params.df <- 
    expand_grid(
        zero.p=c(0.8),
        A.min=c(5)
    )
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )

###################################################
# Generate DAC results for each comparison
###################################################
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
register(MulticoreParam(workers=parsed.args$threads * 2 / 4), default=TRUE)
plan(multisession,      workers=parsed.args$threads * 2 / 4)
# GRanges object with Centro/Telomere regions to filter
data('hg38_cyto') 
# List all separate sample sets + parameters to run multiHiCComapre for
comparisons.df <- 
    ALL_SAMPLE_GROUP_COMPARISONS %>% 
    set_up_sample_comparisons(
        resolutions=parsed.args$resolutions,
        merging='individual'
    ) %>% 
    rename(
        'Sample.Group.P1'=Sample.Group.Numerator,
        'Sample.Group.P2'=Sample.Group.Denominator
    ) %>% 
    run_all_multiHiCCompare(    
        hyper.params.df=hyper.params.df,
        remove.regions=hg38_cyto,
        covariates.df=NULL,
        chromosomes=CHROMOSOMES,
        force_redo=parsed.args$force.redo,
        sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC,
        group1_colname='Sample.Group.P1',
        group2_colname='Sample.Group.P2',
    )

