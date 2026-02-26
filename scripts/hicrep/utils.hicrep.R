###################################################
# Dependencies
###################################################
# library(tidyverse)
# library(magrittr)
library(furrr)
library(glue)

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
    HICREP_RESULTS_DIR %>% 
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
                'SampleID.P1',
                'SampleID.P2'
            )
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
    # Now extract sample metadata from IDs and 
    # format metadata per pair of samples for easy grouping+plotting
    mutate(
        tidy.metadata=
            tidy_pair_metadata(
                sampleID.pairs.df=
                    select(
                        .data=., 
                        all_of(c('SampleID.P1', 'SampleID.P2'))
                    )
            )
    ) %>%
    unnest(tidy.metadata)
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
            # future_pmap(
            pmap(
                .l=.,
                .f=load_hicrep_results,
                .progress=TRUE
            )
    ) %>%
    unnest(hicrep.results)
}

post_proces_hicrep_results <- function(results.df){
    results.df %>% 
    mutate(
        isMerged=
            case_when(
                 str_detect(SampleID.P1, 'Merged') &  str_detect(SampleID.P2, 'Merged') ~ 'Merged vs Merged',
                !str_detect(SampleID.P1, 'Merged') &  str_detect(SampleID.P2, 'Merged') ~ 'Individual vs Merged',
                 str_detect(SampleID.P1, 'Merged') & !str_detect(SampleID.P2, 'Merged') ~ 'Individual vs Merged',
                !str_detect(SampleID.P1, 'Merged') & !str_detect(SampleID.P2, 'Merged') ~ 'Individual vs Individual',
                TRUE ~ NA
            ) %>%
            factor(
                levels=
                    c(
                         'Merged vs Merged',
                         'Individual vs Merged',
                         'Individual vs Individual'
                    )
            )
    ) %>% 
    filter(isMerged != 'Individual vs Merged') %>%
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
    select(-c(filepath))
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

