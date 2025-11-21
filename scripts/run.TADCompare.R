###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/run.TADCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'utils.TADs.R'))
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
# TADCompare parameters
hyper.params.df <- 
    expand_grid(
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )

###################################################
# Generate Consensus TAD calls
###################################################
# List all groups of libraries to call ConsensusTADs() from
sample.groups.df <- 
    tribble(
        ~Sample.Group,
        '16p.iN.WT',
        '16p.iN.DEL',
        '16p.iN.DUP',
        '16p.NSC.WT',
        '16p.NSC.DEL',
        '16p.NSC.DUP'
    ) %>% 
    mutate(
        across(
            starts_with('Sample.Group'),
            ~ str_replace_all(.x, 'All', '.*'),
            .names='{.col}.Pattern'
        )
    ) %>% 
    set_up_sample_groups()
# Call Consensus TADs from individual groups of libraries
sample.groups.df %>%
    run_all_ConsensusTADs(
        hyper.params.df=hyper.params.df,
        # force_redo=TRUE,
        chromosomes=CHROMOSOMES
    )

###################################################
# Generate TADCompare results from merged matrices
###################################################
# List of pairs of merged matrices to compare
merged.matrices.df <- 
    list_mcool_files() %>%
    get_min_resolution_per_matrix() %>% 
    filter(isMerged) %>% 
    mutate(Sample.Group=glue('{Edit}.{Celltype}.{Genotype}')) %>% 
    select(isMerged, resolution, Sample.Group, filepath)
merged.matrices.df
    # nest(data=-c(isMerged, resolution, Sample.Group))
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
    cross_join(tibble(resolution=c(10, 25, 50, 100) * 1e3))
    # cross_join(tibble(resolution=unique(merged.matrices.df$resolution)))
comparisons.df

# merged.matrices.df %>% 
comparisons.df %>% 
    run_all_TADCompare(    
        hyper.params.df=hyper.params.df,
        # force_redo=TRUE,
        chromosomes=CHROMOSOMES
    )

