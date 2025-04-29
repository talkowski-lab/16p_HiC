# Dependencies 
# library(tidyverse)
# library(magrittr)
# library(glue)
# library(purrr)
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
        show_col_types=FALSE
    )
}

load_hiTAD_DIs <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        col_names=
            c(
                'chr',
                'bin.start',
                'bin.end',
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
    } else if (method == 'arrowhead') {
        load_arrowhead_TAD_annotation(filepath)
    } else {
        message(glue('Invalid method {method}'))
    }
}

load_all_TAD_annotations <- function(
    tad_annotations_dir,
    tad.method,
    file_suffix,
    param_names,
    ...){
    tad_annotations_dir %>%
    list.files(
        pattern=glue('*{file_suffix}'),
        full.names=FALSE,
        recursive=TRUE
    ) %>%
    tibble(filename=.) %>% 
    mutate(
        filepath=file.path(tad_annotations_dir, filename),
        filename=str_remove(filename, file_suffix)
    ) %>% 
    separate_wider_delim(
        filename,
        delim=fixed('/'),
        names=param_names
    ) %>%
    mutate(
        TADs=
            pmap(
                .,
                load_TAD_annotation,
                method=tad.method,
                .progress=TRUE
            )
    ) %>%
    select(-c(filepath)) %>%
    unnest(TADs) %>%
    add_column(Method=tad.method)
}

###############
# Compute stuff

###############
# Plot stuff
