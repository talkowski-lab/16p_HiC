###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/compartments/run.compartments.dcHiC.R')
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
})

###################################################
# Set up all hyper-params 
###################################################
options(scipen=999)
# dcHiC hyper-params
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
hyper.params.df <- 
    expand_grid(
        contact.type=c('cis', 'trans'),
        resolution=parsed.args$resolutions
    )

###################################################
# Generate differential compartment results 
###################################################
# List all preprocessed input files per sample group
comparisons.df <- 
    load_dcHiC_sample_groups() %>% 
    setup_dcHiC_group_comparisons(
        comparisons.list=SAMPLE_GROUP_COMPARISONS ,
        cols_to_pair=
            c(
                'resolution',
                'contact.type',
                'Celltype'
            )
    ) %>% 
    mutate(
        output_dir=
            file.path(
                COMPARTMENTS_RESULTS_DIR,
                "method_dcHiC",
                glue('contact.type_{contact.type}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
            ),
        input_filepath=
            file.path(
                output_dir, 
                glue("{Sample.Group.Numerator}-{Sample.Group.Denominator}-dcHiC.input.tsv")
            ),
        script_filepath=
            file.path(
                output_dir, 
                glue("{Sample.Group.Numerator}-{Sample.Group.Denominator}-cmd.sh")
            )
    )
# generate list of commands to run, 1 per input file (pair of matrices)
comparisons.df
comparisons.df %>%
    # add command args
    add_column(
        force_redo=parsed.args$force.redo,
        threads=parsed.args$threads,
        seed=9,
        genome_name="hg38",
        dchic.script.filepath=file.path(SCRIPT_DIR, 'compartments', 'dchicf.r')
    ) %>% 
    # make commands to run dcHiC for each comparison
    pwalk(
        .l=.,
        .f=make_script_file,
        .progress=TRUE
    ) %>% 
    mutate(full.cmd=glue('bash {script_filepath}')) %>% 
    select(full.cmd) %>% 
    # unnest(full.cmd) %>% 
    write_tsv(
        file.path(COMPARTMENTS_DIR, 'run.dcHiC.cmds.txt'),
        col_names=FALSE
    )

