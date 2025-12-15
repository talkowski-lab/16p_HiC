###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/DifferentialContacts/run.multiHiCCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(optparse)
    library(BiocParallel)
})

###################################################
# Set up all comparisons
###################################################
options(scipen=999)
# Individual sample metadata
sample.metadata.df <- 
    load_sample_metadata(filter=TRUE)
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
# GRanges object with Centro/Telomere regions to filter
data('hg38_cyto') 
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$num.cores} core to parallelize'))
register(MulticoreParam(workers=parsed.args$num.cores * 2 / 4), default=TRUE)
plan(multisession,      workers=parsed.args$num.cores * 2 / 4)
# List all separate sample sets + parameters to run multiHiCComapre for
comparisons.df <- 
    tribble(
    ~Sample.Group.Left, ~Sample.Group.Right,
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
    )
comparisons.df
comparisons.df %>% 
    run_all_multiHiCCompare(    
        hyper.params.df=hyper.params.df,
        remove.regions=hg38_cyto,
        covariates.df=NULL,
        # chromosomes=CHROMOSOMES,
        force_redo=parsed.args$force.redo,
        sample_group_priority_fnc=sample_group_priority_fnc_16p
    )

