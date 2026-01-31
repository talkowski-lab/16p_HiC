###################################################
# Dependencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(tictoc)

###################################################
# Load resutls
###################################################
get_smoothing_param <- function(resolution){
    case_when(
        resolution == 1000000 ~  0,
        resolution == 500000  ~  1,
        resolution == 100000  ~  3,
        resolution == 50000   ~  4,
        resolution == 40000   ~  5,
        resolution == 25000   ~ 10,
        resolution == 10000   ~ 20,
        resolution == 5000    ~ 40
    )
}

list_all_hicrep_results_files <- function(samples_to_keep=NULL){
    # List all results files that have been generated
    HICREP_DIR %>% 
    parse_results_filelist(
        suffix='-hicrep.txt',
        filename.column.name='file.pair',
        param_delim='_'
    ) %>%
    # Separate IDs of 2 matrices being compared for each results file
    separate_wider_delim(
        file.pair,
        delim='-',
        names=
            c(
                'SampleInfo.P1.SampleID',
                'SampleInfo.P2.SampleID'
            )
    ) %>% 
    # Extract sample metadata from IDs
    get_info_from_SampleIDs(
        sample_ID_col='SampleInfo.P1.SampleID',
        col_prefix='SampleInfo.P1.',
        keep_id=TRUE
    ) %>% 
    # Repeat for the other sample in each pair
    get_info_from_SampleIDs(
        sample_ID_col='SampleInfo.P2.SampleID',
        col_prefix='SampleInfo.P2.',
        keep_id=TRUE
    ) %>% 
    # keep only specified samples
    {
        if (!is.null(samples_to_keep)) {
            filter(., SampleInfo.P1.SampleID %in% samples_to_keep) %>% 
            filter(SampleInfo.P2.SampleID %in% samples_to_keep)
        } else {
            .
        }
    } %>% 
    # Now format sample metadata per pair for easy grouping+plotting
    pivot_longer(
        starts_with('SampleInfo'),
        names_prefix='SampleInfo.',
        names_to='metadata.field',
        values_to='value'
    ) %>%
    separate_wider_delim(
        metadata.field,
        delim='.',
        names=c('pair.index', 'metadata.field')
    ) %>% 
    pivot_wider(
        names_from='pair.index',
        values_from='value'
    ) %>% 
    mutate(metadata.value=glue('{P1} vs {P2}')) %>% 
    select(-c(P1, P2)) %>%
    pivot_wider(
        names_from=metadata.field,
        values_from=metadata.value
    )
}

load_hicrep_results <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        skip=2,
        progress=FALSE,
        show_col_types=FALSE,
        col_names=c('hicrep.score')
    ) %>% 
    # every line is a score per chromosome in order
    add_column(chr=factor(CHROMOSOMES, levels=CHROMOSOMES))
}

load_all_hicrep_results <- function(samples_to_keep=NULL){
    # Load all files generated from ./scripts/run.hicrep.sh
    list_all_hicrep_results_files(samples_to_keep=samples_to_keep) %>% 
    # Load hicrep scores for each comparison
    mutate(
        hicrep.results=
            future_pmap(
                .l=.,
                .f=load_hicrep_results,
                .progress=TRUE
            )
    ) %>%
    unnest(hicrep.results)
}

post_proces_hicrep_results <- function(results.df){
    results.df %>% 
    rename_with(~ str_replace(.x, '^SampleInfo.', '')) %>% 
    filter(isMerged != 'Individual vs Merged') %>%
    filter(ReadFilter == 'mapq_30 vs mapq_30') %>%
    mutate(
        is.downsampled.fct=
            ifelse(is.downsampled, 'Downsampled', 'Original') %>%
            factor(levels=c('Downsampled', 'Original'))
    ) %>%
    # mutate(
    #     withinGenotype=
    #         case_when(
    #             Genotype %in% c('DEL vs DEL', 'DUP vs DUP', 'WT vs WT') ~ 'Same Genotype',
    #             TRUE ~ 'Different Genotype'
    #         )
    # ) %>% 
    select(-c(ReadFilter, filepath))
}

###################################################
# plotting
###################################################
make_heatmap_plotdf <- function(
    hicrep.results,
    ...){
    hicrep.results %>% 
    group_by(resolution, Sample.ID) %>% 
    summarize(`Mean HiCRep Score`=mean(hicrep.score)) %>%
    ungroup() %>% 
    # Get sample metadata for both samples per pair
    separate_wider_delim(
        Sample.ID,
        delim=' vs ',
        names=c('A.Sample.ID', 'B.Sample.ID')
    ) %>%
    separate_wider_delim(
        c(A.Sample.ID, B.Sample.ID),
        delim=fixed('.'),
        names_sep='.',
        names=c('Edit', 'Genotype', 'SampleNumber', 'Celltype'),
        cols_remove=FALSE
    ) %>%
    rename(
        'A.Sample.ID'=A.Sample.ID.A.Sample.ID,
        'B.Sample.ID'=B.Sample.ID.B.Sample.ID
    ) %>% 
    rename_with(
        .fn=~ str_replace(.x, '(A|B)\\.Sample\\.ID\\.(.*)$', '\\1.\\2'),
        .cols=matches('^(A|B).Sample.ID.')
    )
}

