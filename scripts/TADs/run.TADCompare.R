###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/TADs/run.TADCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
library(hictkR)
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR, 'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADCompare.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

###################################################
# Handle arguments/parameters
###################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# TADCompare hypper-parameters
hyper.params.df <- 
    expand_grid(
        resolutions=parsed.args$resolutions,
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# Load TAD annotations to compare 
TADs.df <- 
    load_all_TAD_results_for_TADCompare()
# List all pairs of matrices to compare
comparisons.list <- 
    ALL_SAMPLE_GROUP_COMPARISONS %>% 
    rename_with(~ str_replace(.x, 'Sample.Group', 'SampleID')) %>% 
    mutate(
        across(
            everything(),
            ~ paste0(.x, '.Merged.Merged')
        )
    )

###################################################
# Generate TADCompare results from merged matrices
###################################################
# List of pairs of merged matrices to compare
comparisons.df <- 
    # List merged matrices
    list_mcool_files() %>%
    filter(isMerged) %>% 
    select(SampleID, filepath) %>%
    # add the TAD boundaries to compare between matrices
    left_join(
        TADs.df,
        by=join_by(SampleID)
    ) %>% 
    # Specify which comparisons to evaluate
    enumerate_pairwise_comparisons(
        sample.group.comparisons=comparisons.list,
        pair_grouping_cols=c('resolution', 'chr', 'TAD.method', 'TAD.params'),
        sampleID_col='SampleID',
        # suffixes=c('.P1', '.P2'),
        suffixes=c('.Numerator', '.Denominator'),
        SampleID.fields=c(NA, 'Celltype', 'Genotype', NA, NA),
        include_merged_col=FALSE
    )
# used by calls to future_pmap() in functions below
# message(glue('using {parsed.args$threads} core to parallelize'))
# plan(multisession, workers=parsed.args$threads)
# Run TADCompare on everything
comparisons.df %>% 
    run_all_TADCompare(    
        hyper.params.df=hyper.params.df,
        TADs.df=TADs.df,
        force_redo=parsed.args$force.redo
    )

