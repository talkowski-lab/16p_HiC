###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/loops/calculate.loop.nesting.levels.R')
BASE_DIR <- here()
# BASE_DIR <- '/data/talkowski/Samples/WAPL_NIPBL/HiC'
suppressPackageStartupMessages({
    library(purrr)
    library(optparse)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})

###################################################
# Set up all comparisons
###################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )

###################################################
# Load loops
###################################################
# Load all loop data to quantify nesting with
loops.df <- 
    # load loop results
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    # anchor.right is the bin start, so change it to bin end to capture that bin in each loop
    mutate(anchor.right=anchor.right + resolution) %>% 
    # nest so each row holds all loops for a given condition (resolution + SampleID + chr)
    nest(
        loops=
            c(
                chr,
                anchor.left,
                anchor.right,
                count,
                length,
                enrichment,
                log10.qval
            )
    )

###################################################
# Generate bed files for loops
###################################################
# Generate .bedpe files for each condition of loops 
loops.df %>% 
    generate_all_loop_bed_files(force_redo=parsed.args$force.redo)

###################################################
# Generate bedtools cmds to calculate loop nesting
###################################################
# Generate bedops commands to calculate how many loops intersect with each genomic bin e.g.

# bin.10      111111111122222222233333333333444444444455555555556666666666
# bin.01      123456789012345678901234567890123456789012345678901234567890
# L1          --------------|==================|--------------------------
# L2          --------------|============|--------------------------------
# L3          --------------|=========|-----------------------------------
# nesting lvl 000000000000003333333333322211111100000000000000000000000000

# First list all bed-like files just listing every bin in the genome at each specified resolution
# This is the reference that we check for overlaps with against our loops 
GENOME_BINS_FILES_DIR %>%
    parse_results_filelist(suffix='bins.tsv') %>%
    select(-c(MatrixID)) %>%
    dplyr::rename('binlist.filepath'=filepath) %>% 
        # {.} -> bin.files.df
    # Now for every set of loops (bed file generated above) generate a 
    # bedtools command to calculate how many loops overlap each genomic bin
    # These commands are all saved to a file and can be run all at once with gnu parallel 
    generate_loop_nesting_calculation_cmds(force_redo=parsed.args$force.redo)
