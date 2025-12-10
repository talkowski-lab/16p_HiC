###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/hicrep/run.hicrep.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'hicrep/utils.hicrep.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(optparse)
})

###################################################
# Arguments
###################################################
# Parse Arguments
parsed.args <- 
    OptionParser() %>%
    add_option(
        c('-f', '--force'),
        action='store_true',
        default=FALSE,
        dest='force.redo'
    ) %>%
    add_option(
        c('-r', '--resolutions'),
        type='character',
        default=paste(c(10, 25, 50, 100) * 1e3, collapse=','),
        dest='resolutions'
    ) %>%
    parse_args(positional_arguments=TRUE) %>%
    {.$options}
parsed.args$resolutions <- 
    parsed.args$resolutions %>%
    str_split(',') %>%
    lapply(as.integer) %>%
    unlist()

###################################################
# Set up param combinations
###################################################
options(scipen=999)
# All combinations of multiHiCCompare hyper-params to test
hyper.params.df <- 
    expand_grid(
        max.window.size=c(100000, 500000, 1000000, 5000000),
        resolution=parsed.args$resolutions
    ) %>%
    mutate(h=get_smoothing_param(resolution)) %>% 
    cross_join(
        tribble(
            ~is.downsampled.flag, ~is.downsampled,
            '',                  FALSE,
            '--bDownSample ',    TRUE
        )
    )

###################################################
# Generate cli commands for all pairs of individual matrices
###################################################
# List all pairs of matrices
list_mcool_files() %>%
    select(SampleID, filepath, isMerged) %>% 
    get_all_row_combinations(
        cols_to_pair=c('isMerged'),
        suffixes=c('.P1', '.P2'),
        keep_self=FALSE
    ) %>%
    # all pairs of matrices X all sets of hyper params
    cross_join(hyper.params.df) %>% 
    mutate(
        output_dir=
            file.path(
                HICREP_DIR,
                glue("resolution_{resolution}"),
                glue("h_{h}"),
                glue("is.downsampled_{is.downsampled}"),
                glue("window.size_{max.window.size}")
            ),
        output_file=
            file.path(
                output_dir,
                glue("{SampleID.P1}-{SampleID.P2}-hicrep.txt")
            ),
        redundant_file=
            file.path(
                output_dir,
                glue("{SampleID.P2}-{SampleID.P1}-hicrep.txt")
            )
    ) %>%
    mutate(
        across(
            c(filepath.P1, filepath.P2, output_dir, output_file, redundant_file),
            ~ str_replace(.x, BASE_DIR, '.')
        )
    ) %>% 
    # ignore results that exists already
    {
        if (parsed.args$force.redo) {
            .
        } else {
            filter(., !(file.exists(output_file) | file.exists(redundant_file)))
        }
    } %>% 
    mutate(
        cmd=glue("hicrep {is.downsampled.flag}--dBPMax {max.window.size} --binSize {resolution} --h {h} {filepath.P1} {filepath.P2} {output_file}")
    ) %>% 
    select(cmd) %>% 
    write_tsv(
        file.path(HICREP_DIR, 'all.hicrep.cmds.txt'),
        col_names=FALSE
    )

