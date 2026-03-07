###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/run.IDR2D.loops.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    source(file.path(BASE_DIR,   'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR,   'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'hicrep/utils.hicrep.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADCompare.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    library(magrittr)
    library(tidyverse)
})

###################################################
# Options
###################################################
options(scipen=999)
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)

###################################################
# HiCRep
###################################################
hicrep.results <- 
    check_cached_results( 
        HICREP_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_hicrep_results,
        samples_to_keep=NULL
    ) %>% 
    post_proces_hicrep_results() %>% 
    standardize_data_cols(skip.isMerged=TRUE)

###################################################
# TADs
###################################################
TADs.df <- 
    check_cached_results(
        results_file=ALL_TAD_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_TAD_results
    ) %>% 
    filter_TAD_results() %T>%
    write_tsv()

###################################################
# Calculate TAD Similarities
###################################################
all.TAD.pairs.df <- 
    TADs.df %>% 
    # select(-c(Celltype, Genotype, SampleID)) %>% 
    nest(boundaries=c(start, end, length)) %>% 
    check_cached_results(
        results_file=ALL_TAD_SIMILARITY_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=calculate_all_boundary_similarities,
        boundaries.df=.,
        only_keep_overlaps=TRUE,
        sample.group.comparisons=
            ALL_SAMPLE_GROUP_COMPARISONS %>%
            rename(
                'Sample.Group.P1'=Sample.Group.Numerator,
                'Sample.Group.P2'=Sample.Group.Denominator,
            ),
        pair_grouping_cols=
            c(
                # 'isMerged',
                # 'weight',
                'resolution',
                'TAD.method',
                'TAD.params',
                'chr'
            )
    ) %>% 
    write_tsv()

###################################################
# TADCompare
###################################################
all.TADCompare.results.df <- 
    check_cached_results(
        results_file=TADCOMPARE_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_TADCompare_results
    ) %>%
    post_process_TADCompare_results() %>% 
    write_tsv()

###################################################
# Loops
###################################################
loops.df <- 
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %T>%
    write_tsv(FILTERED_COOLTOOLS_LOOPS_RESULTS_FILE)

###################################################
# Loop IDR2D 
###################################################
idr2d.results.df <- 
    check_cached_results(
        results_file=ALL_LOOPS_IDR2D_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_IDR2D_results
    ) %>% 
    post_process_IDR2D_results() %>% 
    filter_loop_IDR2D_results() %T>% 
    write_tsv(FILTERED_LOOPS_FILTERED_IDR2D_RESULTS_FILE)


###################################################
# multiHiCCompare
###################################################
differential.contacts.df <- 
    check_cached_results(
        results_file=MULTIHICCOMPARE_RESULTS_FILE,
        force_redo=TRUE,
        results_fnc=load_all_multiHiCCompare_results,
        sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC,
        sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS,
        resolutions=NULL,
        # gw.fdr.threshold=0.1,
        fdr.threshold=0.1,
        nom.threshold=0.05
    ) %>% 
    post_process_multiHiCCompare_results() %>% 
    standardize_data_cols()

###################################################
# Compartments
###################################################
