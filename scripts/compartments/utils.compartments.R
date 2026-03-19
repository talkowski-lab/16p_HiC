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
    ...){
    mutate(
    ) %>% 
    mutate(group=fct_reorder(group, value, .desc=desc))
}

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

