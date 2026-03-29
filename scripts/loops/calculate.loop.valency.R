###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/calculate.loop.valency.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    library(optparse)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
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
# Calculate Loop Valency for each loop anchor
###################################################
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
#  each row is 1 nested set of loop calls per condition + params + context
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
# Also calculate loop valency  i.e. how many loops each anchor is a part of
message('calculating loop valency...')
check_cached_results(
    results_file=ALL_LOOP_VALENCY_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    results_fnc=calculate_all_loop_valency,
    loops.df=loops.df
)

