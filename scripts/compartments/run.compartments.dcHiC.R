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
    library(optparse)
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
# List all pairs of matrices to compare
comparisons.list <- 
    tribble(
        ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        # '16p.iN.DUP',       '16p.iN.DEL', 
        # '16p.NSC.DUP',      '16p.NSC.DEL',
        # '16p.NSC.DUP',      '16p.iN.DUP',
        # '16p.NSC.DEL',      '16p.iN.DEL',
        # '16p.NSC.WT',       '16p.iN.WT',
        # '16p.iN.DUP',       '16p.iN.WT',  
        # '16p.iN.DEL',       '16p.iN.WT',  
        '16p.NSC.DUP',      '16p.NSC.WT',
        '16p.NSC.DEL',      '16p.NSC.WT'
    )
# List all preprocessed input files per sample group
sample.groups.df <- 
comparisons.df <- 
    load_dcHiC_sample_groups() %>% 
    setup_dcHiC_group_comparisons(
        comparisons.list=comparisons.list,
        cols_to_pair=
            c(
                'resolution',
                'contact.type',
                'Celltype'
            )
    ) %>% 
    mutate(
        output.dir=
            file.path(
                COMPARTMENTS_RESULTS_DIR,
                "method_dcHiC",
                glue('contact.type_{contact.type}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
            ),
        input.filepath=
            file.path(
                output.dir, 
                glue("{Sample.Group.Numerator}-{Sample.Group.Denominator}-dcHiC.input.tsv")
            ),
        output.filepath=
            file.path(
                output.dir, 
                glue("{Sample.Group.Numerator}-{Sample.Group.Denominator}-dcHiC.output.tsv")
            )
    )
# generate list of commands to run, 1 per input file (pair of matrices)
comparisons.df
comparisons.df %>%
    # ingore existing results, dont recompute stuff unless explicitly specified
    {
        if (parsed.args$force.redo) {
            .
        } else {
            filter(., !file.exists(output.filepath))
        }
    } %>% 
    # make input files for dcHiC
    rowwise() %>% 
    pwalk(
        .l=.,
        .f=
            function(input.filepath, total.input.file.contents, ...){
                dir.create(
                    dirname(input.filepath),
                    recursive=TRUE,
                    showWarnings=FALSE
                )
                write_tsv(
                    x=total.input.file.contents,
                    file=input.filepath,
                    col_names=FALSE
                )
        },
        .progress=TRUE,
    ) %>% 
    # add command args
    add_column(
        threads=parsed.args$threads,
        seed=9,
        genome_name="hg38",
        dchic.script.filepath=file.path(SCRIPT_DIR, 'compartments', 'dchicf.r')
    ) %>% 
    # make commands to run dcHiC for each comparison
    mutate(
        full.cmd=
            pmap(
                .l=.,
                make_cmd_list,
                .progress=TRUE
            )
    ) %>% 
    select(full.cmd) %>% 
    write_tsv(
        file.path(COMPARTMENTS_DIR, 'run.dcHiC.cmds.txt'),
        col_names=FALSE
    )

