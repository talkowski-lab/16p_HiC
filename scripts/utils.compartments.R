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
run_HiCDOC <- function(
    samples.df,
    resolution,
    chrThreshold,    # Remove chromosomes which are too small to be useful.
    repThreshold,    # Filter sparse replicates to rm uninformative reps with few interactions.
    posThreshold,    # Filter positions (bins) which have too few interactions.
    loessSampleSize,
    ...){
    # row_index=2; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; chrThreshold=tmp$chrThreshold[[row_index]]; repThreshold=tmp$repThreshold[[row_index]]; posThreshold=tmp$posThreshold[[row_index]]; loessSampleSize=tmp$loessSampleSize[[row_index]];
    # data(exampleHiCDOCDataSet)
    row.tmp <- 
    # exampleHiCDOCDataSet %>% 
    HiCDOCDataSetFromCool(
        samples.df$filepath,
        replicates=samples.df$conditions,
        # replicates=samples.df$replicates,
        conditions=samples.df$conditions,
        binSize=resolution
    ) %>% 
    filterSmallChromosomes(threshold=chrThreshold) %>% 
    filterSparseReplicates(threshold=repThreshold) %>% 
    # filterSparseReplicates(threshold=0.99) %>% 
    # filterSparseReplicates(threshold=0.005) %>% 
    # filterWeakPositions(threshold=posThreshold) %>% 
    filterWeakPositions(threshold=0.1) %>% 
    normalizeTechnicalBiases() %>% 
    normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
    detectCompartments()
}

save_HiDOC_results <- function(
    output_dir,
    Sample.Group.Numerator,
    Sample.Group.Denominator,
    ...){
    # parameters(row.tmp)
    compartments(row.tmp)
    differences(row.tmp)
    concordances(row.tmp)
}

run_all_HiCDOC <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    ...){
    # sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL; force_redo=FALSE;
    comparisons.df %>% 
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
                Sample.Group.P1.Priority < Sample.Group.P2.Priority ~ Sample.Group.P1,
                Sample.Group.P1.Priority > Sample.Group.P2.Priority ~ Sample.Group.P2,
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
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        output_dir=
            file.path(
                COMPARTMENTS_DIR,
                'results',
                'method_HiCDOC',
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('chrThreshold_{chrThreshold}'),
                glue('repThreshold_{repThreshold}'),
                glue('posThreshold_{posThreshold}'),
                glue('loessSampleSize_{loessSampleSize}')
            ),
    ) %>% 
    arrange(resolution) %>% 
        {.} -> tmp
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
