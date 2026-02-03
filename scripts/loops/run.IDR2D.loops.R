###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/run.IDR2D.loops.R')
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
# All combinations of multiHiCCompare hyper-params to test
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
hyper.params.df <- 
    tribble(
        ~metric_colname, ~value_transformation,
        'enrichment',    'identity',  # high enrichment => most important loops 
        'log10.qval',    'identity'   # high log10.qval => most important loops
    ) %>% 
    cross_join(tibble(ambiguity_resolution_method=c("overlap", "midpoint", "value")))
q.thresh <- 0.1 

###################################################
# Generate DAC results for each comparison
###################################################
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
#  Prepare loops for comparison between conditions 
#  each row is 1 nested set of loop calls per condition + context
nested.loops.df <- 
    # load loop results
    check_cached_results(
        results_file=COOLTOOLS_LOOPS_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    post_process_cooltools_dots_results() %>% 
    # standardize_data_cols() %>% 
    # filter relevant results
    filter(weight == 'balanced') %>% 
    filter(kernel == 'donut') %>% 
    filter(log10.qval > -log10(q.thresh)) %>% 
    # prep columns for input to IDR2D
    mutate(resolution=scale_numbers(resolution, force_numeric=TRUE)) %>% 
    select(-c(Edit, Celltype, Genotype)) %>% 
    rename(
        'start.A'=anchor.left,
        'start.B'=anchor.right
    ) %>% 
    mutate(
        end.A=start.A + resolution,
        end.B=start.B + resolution,
        chr.A=chr,
        chr.B=chr
    ) %>% 
    # nest so one set of loop calls per row (SampleID + res + chr + weight)
    nest(
        loops=
            c(
                chr.A, start.A, end.A,
                chr.B, start.B, end.B, 
                count, length, enrichment, log10.qval
            )
    )
# run IDR2D on all comparisons of sample groups + param sets
nested.loops.df %>% 
    run_all_IDR2D_analysis(
        hyper.params.df=hyper.params.df,
        force.redo=parsed.args$force.redo,
        sample.group.comparisons=
            ALL_SAMPLE_GROUP_COMPARISONS %>% 
            rename(
                'SampleID.P1'=Sample.Group.Numerator,
                'SampleID.P2'=Sample.Group.Denominator
            ),
        pair_grouping_cols=
            c(
                'weight',
                'resolution',
                'kernel',
                'chr'
            ),
        SampleID.fields=
            c(
                NA,
                'Celltype',
                'Genotype'
            )
    )
