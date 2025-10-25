###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(furrr)
library(HiCDOC)

###################################################
# Generate Compartment Results
###################################################
# set_up_sample_comparisons <- function(comparison.groups){
#     # get info + filepaths for all contact matrices
#     load_mcool_files(
#         return_metadata_only=TRUE,
#         keep_metadata_columns=FALSE
#     ) %>% 
#     get_min_resolution_per_matrix() %>% 
#     distinct() %>% 
#     # Now group samples by condition, 
#     filter(!isMerged) %>% 
#     nest(samples.df=-c(isMerged)) %>% 
#     cross_join(comparison.groups) %>%
#     # subset relevant samples for each comparison
#     rowwise() %>% 
#     mutate(
#         samples.df=
#             samples.df %>%
#             mutate(
#                 Sample.Group=
#                     case_when(
#                         str_detect(SampleID, Sample.Group.P1.Pattern) ~ Sample.Group.P1,
#                         str_detect(SampleID, Sample.Group.P2.Pattern) ~ Sample.Group.P2,
#                         TRUE ~ NA
#                     )
#             ) %>%
#             filter(!is.na(Sample.Group)) %>% 
#             list()
#     ) %>%
#     # minimum and max resoltion of all individual matrices per comparison
#     mutate(
#         resolution.min=min(samples.df$resolution),
#         resolution.max=max(samples.df$resolution)
#     ) %>% 
#     # Now list every comparison at every resolution that is either a min or max for 1 comparison
#     ungroup() %>%
#     mutate(resolution=list(unique(c(resolution.min, resolution.max)))) %>%
#     unnest(resolution) %>% 
#     mutate(
#         resolution.type=
#             case_when(
#                 resolution == resolution.max ~ 'max',
#                 resolution == resolution.min ~ 'min',
#                 TRUE                         ~ NA
#             )
#     ) %>% 
#     select(-c(resolution.min, resolution.max))
# }
run_HiCDOC <- function(
    samples.df,
    sample_group_priority_fnc,
    filepath,
    resolution,
    replicates,
    conditions,
    # range1,
    # range2,
    effect.col='conditions',
    replicate.col='replicates'
    ...){
    # data(exampleHiCDOCDataSet)
    HiCDOCDataSetFromCool(
        filepath,
        replicates=replicate,
        conditions=condition,
        binSize=resolution
    ) %>% 
    filterSmallChromosomes(threshold=100) %>% 
    filterSparseReplicates(threshold=0.3)
    filterWeakPositions(threshold=1)
    normalizeTechnicalBiases()
    normalizeDistanceEffect(loessSampleSize=20000)
    detectCompartments()
}

run_all_HiCDOC <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    covariates.df=NULL,
    chromosomes=CHROMOSOMES,
    ...){
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(tibble(chr=chromosomes)) %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    mutate(
        across(
            matches('Sample.Group.(P1|P2)'),
            ~ sample_group_priority_fnc(.x),
            .names='{.col}.Priority'
        )
    ) %>%
    mutate(
        Sample.Group.Numerator=
            case_when(
                Sample.Group.P1.Priority > Sample.Group.P2.Priority ~ Sample.Group.P1,
                Sample.Group.P1.Priority < Sample.Group.P2.Priority ~ Sample.Group.P2,
                TRUE ~ NA
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == Sample.Group.P1 ~ Sample.Group.P2,
                Sample.Group.Numerator == Sample.Group.P2 ~ Sample.Group.P1,
                TRUE ~ NA
            )
    ) %>% 
    select(-c(matches('Sample.Group.(P1|P2)'))) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                COMPARTMENTS_DIR,
                'results',
                'method_HiCDOC',
                glue('merged_{isMerged}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-compartments.tsv')
            )
    ) %>% 
    arrange(resolution) %>% 
        # {.} -> tmp
    # nibpl.tmp <- tmp %>% filter(Sample.Group.Numerator == "NIPBL.iN.DEL", Sample.Group.Denominator == "NIPBL.iN.WT") %>% select(samples.df, Sample.Group.Numerator, Sample.Group.Denominator, chr)
    # nibpl.tmp$samples.df[[1]]
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_HiCDOC,
                    sample_group_priority_fnc=sample_group_priority_fnc,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load Compartment Results
###################################################
