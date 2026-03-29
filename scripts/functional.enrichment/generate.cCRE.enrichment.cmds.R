###################################################
# Dependencies
###################################################
library(here)
here::i_am('notebooks/TADs.Rmd')
BASE_DIR <- here()
# BASE_DIR <- '/data/talkowski/Samples/WAPL_NIPBL/HiC'
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'functional.enrichment/utils.enrichment.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    library(tidyverse)
    library(magrittr)
})
options(scipen=999)
options(future.globals.maxSize=2.5 * (1024 ** 3))
# plan(sequential)
plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)

###################################################
# CLI Args
###################################################
options(scipen=999)
options(show_col_types=FALSE)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
FORCE_REDO <- parsed.args$force.redo

###################################################
# Generate cCRE BED files BED files
###################################################
# First separate ENCODE cCREs into separate bed files per type
load_encode_ccres() %>% 
    split_cCRE_annotations()

###################################################
# Generate TADs BED files
###################################################
# Load all TAD annotations
message('Loading TADs...')
TADs.df <- 
    ALL_TAD_RESULTS_FILE %>%
    read_tsv() %>% 
    filter(TAD.length > 0) %>% 
    mutate(
        param_dir=
            file.path(
                glue('method_{method}'),
                glue('TAD.params_{TAD.params}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('region_{chr}')
            )
    )
message('Generating TAD BED files...')
# Make bed files for the span of each TAD in a nested directory structure for bedtools queries
TADs.df %>% 
    format_TADs_for_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Make bed files for just the TAD anchors in a nested directory structure for bedtools queries
TADs.df %>% 
    format_TAD_anchors_for_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)

rm(TADs.df)
message('Finished generating TAD BED files')
###################################################
# Generate Loop BED files
###################################################
# Load all loop annotations
message('Loading loops...')
loops.df <- 
    ALL_COOLTOOLS_LOOPS_RESULTS_FILE %>%
    read_tsv() %>% 
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    mutate(
        param_dir=
            file.path(
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('region_{chr}')
            )
    )
message('Generating all loop BED files...')
# Make bed files for the span of each loop in a nested directory structure for bedtools queries
loops.df %>% 
    make_loop_span_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Make bed files for just the loop anchors in a nested directory structure for bedtools queries
loops.df %>% 
    make_loop_bedpe_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Make bed files for just the loop anchors in a nested directory structure for bedtools queries
# We already generated this data with the loop valency analysis, so we can use that and convert it to structured bed files.
ALL_LOOP_VALENCY_RESULTS_FILE %>%
    read_tsv() %>% 
    post_process_loop_valency_results() %>% 
    make_loop_anchor_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Also make bed files for the loop nesting regions, where every row represents a contiguous set of bins that are overlapped by the exact same set of loops.
ALL_LOOP_NESTING_RESULTS_FILE %>%
    read_tsv() %>% 
    make_loop_nesting_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)

rm(loops.df)
message('Finished generating all loop BED files')
###################################################
# Generate Compartments BED files
###################################################
# Load compartment data
# compartments.df <- 
# Compartment interiors
# Compartment boundaries

# rm(compartments.df)
###################################################
# Generate DIRs BED files
###################################################
# Load all DIRs
message('Loading DIRs...')
DIRs.df <- 
    ALL_MULTIHICCOMPARE_RESULTS_FILE %>%
    read_tsv() %>% 
    post_process_multiHiCCompare_results() %>% 
    select(-c(bin.pair.idx, distance.bp)) %>% 
    mutate(
        param_dir=
            file.path(
                glue('merged_{merged}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
            )
    )
message('Generating DIR BED files...')
# Make bed files for the region between DIR anchors in a nested directory structure for bedtools queries. Unlike with TADs or loops, I do not necessairily expect the "inside" of a DIR to show signal unless it overlps with some other features.
DIRs.df %>% 
    make_DIR_span_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Make bedpe files, which store both anchors on the same line, to preserve the association
DIRs.df %>% 
    make_DIR_bedpe_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)
# Make bed files for just the DIR anchors in a nested directory structure for bedtools queries
DIRs.df %>% 
    make_DIR_anchor_bed_files() %>% 
    write_all_bed_files(force_redo=FORCE_REDO)

rm(DIRs.df)
message('Finished generating DIR BED files')
###################################################
# Generate bedtools intersect cmds for cCRE enrichment
###################################################
# First list all bed files for all annotated feautres
list(
    bed_dir=
        c(
            TAD_BED_FILES_DIR,
            LOOP_BED_FILES_DIR,
            # COMPARTMENT_BED_FILES_DIR,
            MULTIHICCOMPARE_BED_FILES_DIR
        )
) %>% 
# List all bed files in all feature directories
pmap(
    .f=list_bed_files,
    suffix='-.*.bed(pe)?'
) %>%
bind_rows() %>% 
filter(
    region != 'DIR.spans',
    file.type == 'bed'
) %>% 
# Next generate `bedtools intersect` commands to calculate how many cCRES of each type overlap with the feature region or just at the anchors/sites/boundaries.
make_bedtools_cmds(
    force_redo=FORCE_REDO,
    bed_cmds_file=ALL_CCRE_ENRICHMENT_CMDS_FILE
)
message(glue('Generated file with all cCRE enrichment commands: {ALL_CCRE_ENRICHMENT_CMDS_FILE}'))
# Now you can run the commands in the file with `gnu parallel` to generate all the enrichment results files (1 per line)

