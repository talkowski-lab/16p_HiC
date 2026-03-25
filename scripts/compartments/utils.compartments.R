###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(furrr)

###################################################
# Generate cooltools results
###################################################
quantize_track <- function(
    values,
    thresholds,
    desc=FALSE,
    ...){
    values %>% 
    tibble(value=.) %>%
    # thresholds must be a named list
    mutate(
        group=
            cut(
                value,
                breaks=c(-Inf, thresholds, Inf),
                labels=
                    ifelse(
                        is.null(names(thresholds)),
                        paste('<', thresholds),
                        names(thresholds)
                    )
            )
    ) %>% 
    mutate(group=fct_reorder(group, value, .desc=desc)) %>%
    select(group)
}

load_cooltools_compartment_results <- function(filepath){
    filepath %>%
    read_tsv(show_col_types=FALSE)
}

load_all_cooltools_compartment_results <- function(){
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix='-.cis.vecs.tsv'
    ) %>%  
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
    results.df %>%
    quantize_track_values(
        column_to_quantize='eigs',
        n_quantiles=5,
        labels=c('A', 'B'),
        thresholds=c(0),
        thresh.min=-Inf,
        thresh.max=Inf
    )
}

