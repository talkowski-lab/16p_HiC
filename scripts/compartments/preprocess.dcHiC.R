###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/compartments/preprocess.dcHiC.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
    library(optparse)
})

###################################################
# Set up all hyper-params 
###################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# All combinations of multiHiCCompare hyper-params to test
hyper.params.df <- 
    expand_grid(
        contact.type=c('cis', 'trans'),
        resolution=parsed.args$resolutions
    )

###################################################
# Generate input files for dcHiC
###################################################
PREPROCESS_SCRIPT_FILEPATH <- file.path(SCRIPT_DIR, 'compartments', 'preprocess.dcHiC.py')
list_mcool_files() %>% 
    filter(!isMerged) %>% 
    select(SampleID, filepath) %>% 
    cross_join(hyper.params.df) %>% 
    mutate(
        hyper.param.path=
            file.path(
                glue("contact.type_{contact.type}"),
                glue("resolution_{resolution}")
            ),
        preprocess.output.dir=
            file.path(
                COMPARTMENTS_PREPROCESSED_DIR,
                hyper.param.path
            ),
        preprocess.cmd=glue("mkdir -p {preprocess.output.dir}; python {PREPROCESS_SCRIPT_FILEPATH} -genomeFile {CHROMOSOME_SIZES_FILE} -res {resolution} -prefix {SampleID} -input cool -output {preprocess.output.dir} -file {filepath}"), #%>% str_replace_all(BASE_DIR, '.'),
        preprocessed.matrix.filepath=
            file.path(
                COMPARTMENTS_PREPROCESSED_DIR,
                hyper.param.path,
                glue("{SampleID}.matrix")
            ),
        preprocessed.region.filepath=
            file.path(
                COMPARTMENTS_PREPROCESSED_DIR,
                hyper.param.path,
                glue("{SampleID}.abs.bed")
            )
    ) %>% 
    {
        if (parsed.args$force.redo) {
            .
        } else {
            filter(., !file.exists(preprocessed.matrix.filepath))
        }
    } %>% 
    select(preprocess.cmd) %>% 
    # head(1) %>% t()
    write_tsv(
        file.path(COMPARTMENTS_DIR, 'preprocess.dcHiC.cmds.txt'),
        col_names=FALSE
    )
