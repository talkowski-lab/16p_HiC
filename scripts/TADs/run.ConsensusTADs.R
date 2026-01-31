###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/TADs/run.TADCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(hictkR)
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR, 'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADCompare.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
})

###################################################
# Handle arguments/parameters
###################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# TADCompare parameters
hyper.params.df <- 
    expand_grid(
        resolutions=parsed.args$resolutions,
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )

###################################################
# Generate ConsensusTAD results
###################################################
comparisons.df <- 
    ALL_SAMPLE_GROUPS %>% 
    set_up_sample_groups(
        resolutions=parsed.args$resolutions,
        use_merged=FALSE
    )
# used by calls to future_pmap() in functions below
if (parsed.args$threads > 1) {
    message(glue('using {parsed.args$threads} core to parallelize'))
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}
# Generate results
# comparisons.df
comparisons.df %>% 
    run_all_ConsensusTADs(
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo,
        chromosomes=CHROMOSOMES  # all chromosomes
    )

