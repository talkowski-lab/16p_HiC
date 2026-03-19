###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(furrr)
library(HiCDOC)

###################################################
# Generate HiCDOC Results
###################################################
run_HiCDOC <- function(
    samples.df,
    resolution,
    chrThreshold,    # Remove chromosomes which are too small to be useful.
    repThreshold,    # Filter sparse replicates to rm uninformative reps with few interactions.
    posThreshold,    # Filter positions (bins) which have too few interactions.
    loessSampleSize,
    ...){
    # row_index=1; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; chrThreshold=tmp$chrThreshold[[row_index]]; repThreshold=tmp$repThreshold[[row_index]]; posThreshold=tmp$posThreshold[[row_index]]; loessSampleSize=tmp$loessSampleSize[[row_index]];
    tmp1 <- 
        HiCDOCDataSetFromCool(
            paths=samples.df$filepath,
            replicates=samples.df$replicates,
            conditions=samples.df$conditions,
            binSize=resolution
        )
    tmp2 <- 
        tmp1 %>% 
        filterSmallChromosomes(threshold=chrThreshold) %>% 
        # filterSparseReplicates(threshold=repThreshold) %>% 
        filterSparseReplicates(threshold=0.005) %>% 
        # filterWeakPositions(threshold=posThreshold) %>% 
        filterWeakPositions(threshold=0.01) %>% 
        normalizeTechnicalBiases() %>% 
        normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
        detectCompartments()
    names(tmp2)
    str(tmp2)
    compartments(tmp2)
    differences(tmp2)
    concordances(tmp2)
    hicdoc.results.obj <- 
        HiCDOCDataSetFromCool(
            samples.df$filepath,
            replicates=samples.df$replicates,
            conditions=samples.df$conditions,
            binSize=resolution
        ) %>% 
        filterSmallChromosomes(threshold=chrThreshold) %>% 
        filterSparseReplicates(threshold=repThreshold) %>% 
        # filterSparseReplicates(threshold=0.99) %>% 
        # filterSparseReplicates(threshold=0.005) %>% 
        filterWeakPositions(threshold=posThreshold) %>% 
        # filterWeakPositions(threshold=0.1) %>% 
        normalizeTechnicalBiases() %>% 
        normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
        detectCompartments()
    compartments(hicdoc.results.obj)
    differences(hicdoc.results.obj)
    concordances(hicdoc.results.obj)
}

run_all_HiCDOC <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    ...){
    sample_group_priority_fnc=sample_group_priority_fnc_Cohesin; force_redo=FALSE;
    # sample_group_priority_fnc=sample_group_priority_fnc_16p; force_redo=FALSE;
    comparisons.df %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    set_foldchange_direction_as_factor(
        sample_group_priority_fnc=sample_group_priority_fnc,
        priority_col_pattern=starts_with('Sample.Group.'),
        label_cols=c('Sample.Group.Left', 'Sample.Group.Right')
    ) %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        output_dir=
            file.path(
                COMPARTMENTS_RESULTS_DIR,
                'method_HiCDOC',
                glue('resolution_{scale_numbers(resolution)}'),
                glue('chrThreshold_{chrThreshold}'),
                glue('repThreshold_{repThreshold}'),
                glue('posThreshold_{posThreshold}'),
                glue('loessSampleSize_{loessSampleSize}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-HiCDOC.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !(file.exists(results_file)))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            dplyr::count(
                .,
                resolution,
                chrThreshold,
                repThreshold,
                posThreshold,
                loessSampleSize,
                Sample.Group.Numerator,
                Sample.Group.Denominator
            ) %>%
            rename_with(~ str_remove(., 'Sample.Group.'))
        )
    } %>%
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
                    ...  # passed from the call run_HiCDOC()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Generate cooltools results
###################################################
load_cooltools_compartment_results <- function(filepath){
    filepath %>%
    read_tsv(show_col_types=FALSE)
}

list_all_cooltools_compartment_results <- function(){
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        suffix='-cis.vecs.tsv',
        filename.column.name='MatrixID'
    ) # %>% 
    # get_info_from_MatrixIDs(keep_id=FALSE)
}

load_all_cooltools_compartment_results <- function(){
    list_all_cooltools_compartment_results() %>%
    mutate(
        compartments=
            future_pmap(
                 .l=.,
                 .f=load_cooltools_compartment_results,
                 .progress=TRUE
            )
    ) %>%
    dplyr::rename('chr'=chrom) %>%
    select(-c(filepath))

}

post_process_cooltools_compartment_results <- function(results.df){
}

