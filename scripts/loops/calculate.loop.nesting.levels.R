###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/calculate.loop.nesting.levels.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    library(optparse)
    source(file.path(BASE_DIR,   'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR,   'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'loops', 'utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})

###################################################
# Set up all comparisons
###################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )

###################################################
# Generate DAC results for each comparison
###################################################
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
#  Prepare loops for comparison between conditions 
#  each row is 1 nested set of loop calls per condition + context
loops.df <- 
    # load loop results
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    # standardize_data_cols() %>% 
    nest(
        loops=
            c(
                anchor.left,
                anchor.right,
                count,
                length,
                enrichment,
                log10.qval
            )
    )
# For each bin calculate how many loops overlap it i.e. how "nested" that bin is
check_cached_results(
    results_file=ALL_LOOP_NESTING_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    results_fnc=calculate_all_loop_nesting,
    loops.df=loops.df
)

