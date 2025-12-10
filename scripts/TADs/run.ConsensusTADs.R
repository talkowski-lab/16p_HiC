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
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADCompare.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(furrr)
    library(optparse)
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
# Generate ConsensusTAD results
###################################################
# TADCompare parameters
hyper.params.df <- 
    expand_grid(
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# used by calls to future_pmap() in functions below
message(glue('using {parsed.args$num.cores} core to parallelize'))
plan(multisession, workers=parsed.args$num.cores)
# Set up results to generate
tribble(
    ~Sample.Group,
    # '16p.iN.DUP',
    # '16p.iN.DEL',
    # '16p.iN.WT',
    '16p.NSC.DUP',
    '16p.NSC.DEL',
    '16p.NSC.WT'
) %>% 
set_up_sample_groups(
    resolutions=parsed.args$resolutions,
    use_merged=FALSE
) %>% 
run_all_ConsensusTADs(
    hyper.params.df=hyper.params.df,
    force_redo=parsed.args$force.redo,
    chromosomes=CHROMOSOMES  # all chromosomes
)

