library(tidyverse)
library(stringi)
# library(HiContacts)
library(glue)
library(furrr)
library(ggplot2)
library(viridis)
library(ggh4x)
# library(ggtext)
library(cowplot)
###############
# Load TAD annotation data
load_arrowhead_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        skip=2,
        show_col_types=FALSE,
        col_select=c(1,2,3,12,13,14),
        # col_select=c(1,2,3,4,5,6,12,13,14),
        col_names=
            c(
                'chr',
                'start',
                'end',
                'B.chr',
                'B.start',
                'B.end',
                'name',
                'score',
                'strand1',
                'strand2',
                'color',
                'corner.score',
                'uVarScore',
                'lVarScore',
                'upSign',
                'loSign'
            )
    )
}

load_hiTAD_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'TAD.start',
                'TAD.end'
            )
    )
}

load_cooltools_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        progress=FALSE,
        show_col_types=FALSE
    ) %>%
    rename(
        'chr'=chrom,
        'TAD.start'=start,
        'TAD.end'=end
    ) %>% 
    select(-c(region)) %>% 
    pivot_longer(
        ends_with('00'),
        names_to='TAD.stat',
        values_to='value'
    ) %>%
    mutate(TAD.stat=stri_reverse(TAD.stat)) %>% 
    separate_wider_delim(
        TAD.stat,
        delim='_',
        names=c(
            'window.size',
            'stat'
        ),
        too_many='merge'
    ) %>%
    mutate(
        across(
            c(stat, window.size),
            stri_reverse
        ),
        window.size=as.integer(window.size)
    )
}

load_hiTAD_DIs <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'start',
                'end',
                'DI'
            )
    )
}

load_TAD_annotation <- function(
    filepath,
    method,
    ...){
    if (method == 'hiTAD') {
        load_hiTAD_TAD_annotation(filepath)
    } else if (method == 'hiTAD-DIs') {
        load_hiTAD_DIs(filepath)
    } else if (method == 'cooltools') {
        load_cooltools_TAD_annotation(filepath)
    } else if (method == 'arrowhead') {
        load_arrowhead_TAD_annotation(filepath)
    } else {
        message(glue('Invalid method {method}'))
    }
}

load_all_TAD_annotations <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-TAD.tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    process_matrix_name() %>% 
    mutate(
        TADs=
            pmap(
                .,
                load_TAD_annotation,
                .progress=TRUE
            )
    ) %>%
    select(-c(filepath)) %>% 
    unnest(TADs)
}
###############
# Compute stuff

###############
# Plot stuff
