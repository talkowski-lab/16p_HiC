###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/TADs/run.TADCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

###################################################
# Arguments
###################################################
# Parse Arguments
parsed.args <- 
    OptionParser() %>%
    add_option(
        c('-t', '--threads'),
        type='integer',
        default=length(availableWorkers()),
        dest='num.cores'
    ) %>% 
    add_option(
        c('-f', '--force'),
        action='store_true',
        default=FALSE,
        dest='force.redo'
    ) %>%
    add_option(
        c('-r', '--resolutions'),
        type='character',
        default=paste(c(10, 25, 50, 100) * 1e3, collapse=','),
        dest='resolutions'
    ) %>%
    parse_args(positional_arguments=TRUE) %>%
    {.$options}
parsed.args$resolutions <- 
    parsed.args$resolutions %>%
    str_split(',') %>%
    lapply(as.integer) %>%
    unlist()

###################################################
# Generate TADCompare results from merged matrices
###################################################
# TADCompare hypper-parameters
hyper.params.df <- 
    expand_grid(
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# List of pairs of merged matrices to compare
merged.matrices.df <- 
    list_mcool_files() %>%
    get_min_resolution_per_matrix() %>% 
    filter(isMerged) %>% 
    mutate(Sample.Group=glue('{Edit}.{Celltype}.{Genotype}')) %>% 
    select(isMerged, resolution, Sample.Group, filepath)
# get all relevant pairs of merged matrices
comparisons.df <- 
    tribble(
        ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        '16p.iN.DEL',            '16p.iN.WT',  
        '16p.iN.DUP',            '16p.iN.WT',  
        # '16p.iN.DUP',            '16p.iN.DEL', 

        '16p.NSC.DEL',           '16p.NSC.WT',
        '16p.NSC.DUP',           '16p.NSC.WT',
        # '16p.NSC.DUP',           '16p.NSC.DEL',

        # '16p.iN.DEL',            '16p.NSC.DEL',
        # '16p.iN.DUP',            '16p.NSC.DUP',
        '16p.iN.WT',             '16p.NSC.WT'
    ) %>% 
    left_join(
        merged.matrices.df %>% select(-c(resolution)),
        by=join_by(Sample.Group.Numerator == Sample.Group)
    ) %>% 
    left_join(
        merged.matrices.df %>% select(-c(resolution)),
        suffix=c('.Numerator', '.Denominator'),
        by=join_by(isMerged, Sample.Group.Denominator == Sample.Group)
    ) %>% 
    cross_join(tibble(resolution=parsed.args$resolutions))
# Run TADCompare on everything
message(glue('using {parsed.args$num.cores} core to parallelize'))
plan(multisession, workers=parsed.args$num.cores)
comparisons.df %>% 
    run_all_TADCompare(    
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo,
        chromosomes=CHROMOSOMES
    )

