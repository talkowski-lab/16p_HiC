###################################################
# Dependencies 
###################################################
library(here)
here::i_am('scripts/compartments/run.compartment.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts',      'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    # source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'compartments', 'utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(future)
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
# Hyper-params
###################################################
# List hyper parameter combinations to use
hyper.params.df <- 
    expand_grid(
        chrThreshold=c(100),
        repThreshold=c(0.05),
        posThreshold=c(1),
        loessSampleSize=c(20000)
    )

###################################################
# Load data + Set up comparisons to compute 
###################################################
individual.matrix.files <- 
    list_mcool_files(pattern='.hg38.mapq_30.1000.mcool') %>%
    filter(!isMerged) %>% 
    mutate(condition=glue('{Edit}.{Celltype}.{Genotype}')) %>% 
    mutate(replicate=glue('{CloneID}.{Celltype}.{Genotype}')) %>%
    cross_join(tibble(resolution=parsed.args$resolution)) %>% 
    nest(data=-c(condition, resolution))

###################################################
# Generate Compartments
###################################################
comparisons.df <- 
    tribble(
        ~Sample.Group.Left, ~Sample.Group.Right,
    ) %>% 
    set_up_sample_comparisons(
        resolutions=parsed.args$resolutions,
        merging='individual'
    )
    select(-c(ends_with('.Pattern'))) %>%
    rowwise() %>%
    mutate(
        samples.df=
            samples.df %>%
            dplyr::rename('conditions'=Sample.Group) %>%
            mutate(replicates=factor(CloneID)) %>%
            select(-c(resolution)) %>% 
            list()
    ) %>% 
    ungroup()
# parallelizing params
plan(multisession, workers=parsed.args$threads)
# generate compartments
comparisons.df %>% 
    run_all_HiCDOC(    
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo,
        sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL
    )
