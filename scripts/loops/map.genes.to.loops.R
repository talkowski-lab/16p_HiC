###################################################
# Dependencies
###################################################
library(here)
here::i_am('scripts/loops/map.genes.to.loops.R')
BASE_DIR <- here()
source(file.path(BASE_DIR, 'scripts/locations.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(BASE_DIR, 'scripts/constants.R'))
source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
library(tidyverse)
library(magrittr)

###################################################
# Load gene coords + metadata
###################################################
genes.df <- 
    load_gene_expression_data() %>%
    filter(feature == "protein_coding") %>% 
    select(-c(feature, SampleID, TPM.mean, TPM.sd)) %>% 
    distinct()

###################################################
# Load loop reproducibility results
###################################################
# Load IDR2D results for all comparisons of sample groups
idr2d.results.df <- 
    check_cached_results(
        results_file=LOOPS_IDR2D_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_IDR2D_results
    ) %>% 
    # only use reuslts based on the donut kernel
    filter(kernel == 'donut') %>% 
    # only keep results using 'balanced' matrices as input
    filter(weight == 'balanced') %>% 
    # IDR2D tie-breaking method
    filter(resolve.method == 'value') %>% 
    # rank loops by log10(qvalue) for IDR
    filter(metric == 'log10.qval') %>% 
    mutate(
        is.loop.shared=
            case_when(
                loop.type == 'P1.only' ~ loop.type,
                loop.type == 'P2.only' ~ loop.type,
                IDR < 0.1              ~ 'IDR < 0.1',
                IDR <= 1               ~ 'Irreproducible',
                TRUE                   ~ NA
            )
    ) %>%
    standardize_data_cols()

###################################################
# mapp genes to which loops they overlap
###################################################
genes.df %>% 
    inner_join(
        idr2d.results.df,
        by=
            join_by(
                chr,
                # any gene that overalps the loop
                # overlaps(start, end, anchor.left, anchor.right)
                # any gene that is strictly within the loop
                within(start, end, anchor.left, anchor.right)
            )
    ) %>% 
    select(
        # weight,
        resolution,
        # kernel,
        # metric,
        # resolve.method,
        # comparison,
        SampleID.P1,
        SampleID.P2,
        anchor.left,
        anchor.right,
        diff.value,
        diff.rank,
        IDR,
        # loop.type,
        is.loop.shared,
        EnsemblID,
        symbol,
        chr,
        start,
        end
    ) %T>%
    write_tsv(LOOPS_IDR2D_GENE_MAPPING_FILE)
