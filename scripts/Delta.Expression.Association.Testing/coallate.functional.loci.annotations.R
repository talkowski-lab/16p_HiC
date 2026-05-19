################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/Delta.Expression.Association.Testing/make.clean.functional.locus.gene.mappings.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'FGE.Association.Testing/utils.enrichment.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.association.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
})

################################################################################
# Params to specify
################################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
# parsed.args$threads <- 1
if (parsed.args$threads > 1){
    # plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)
    options(future.globals.maxSize=1.23 * 1024**3)
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}

################################################################################
# Load Association data i.e. functional associations between loci and Genes
################################################################################
# ABC-based enhancers-promoter pairs
enhancers.df <- 
    bind_rows(
        # load_internal_ABC_enhancers() %>% 
        #     add_column(association.source='Internal'),
        load_nasser_ABC_enhancers() %>%
            select(-c(CellType)) %>% 
            add_column(association.source='Nasser_2021')
    ) %>% 
    add_column(association.type='ABC.enhancer')
# eQTLs linking genes to loci with specific variants
# eQTLs.df <- 
#     load_all_clean_eQTLs(parsed.args$fore.redo) %>% 
#     dplyr::rename('association.subtype'=) %>% 
#     add_column(association.source='EBI') %>% 
#     add_column(association.type='eQTLs')
# # TFs linking genes to the locations of TFs genes that target them
# TFs.df <- 
#     tibble()
# TF binding sites, linking genes to bindings sites of TFs that target them
TFBS.df <- 
    load_TF_binding_sites() %>% 
    add_count(TF.Symbol) %>% filter(n > 2000) %>% select(-c(n)) %>% 
    select(
        -c(
            ID,
            antibody.set,
            cell.set,
            exp.set,
            treatment.set,
            `peak-caller.set`,
            `peak-caller.count`,
            `peak-caller.list`,
            # peak.count,
            # summit,
            TF.ClassID,

        )
    ) %>% 
    mutate(Target.Gene.Symbol=TF.Symbol) %>% 
    dplyr::rename('association.subtype'=TF.Symbol) %>% 
    add_column(association.source='GTRD') %>% 
    add_column(association.type='TFBS')

################################################################################
# Combine all functional elements and add gene metadata for downstream analysis
################################################################################
# positions of genes, used for associations with condition-specific HiFs
gene.positions.df <- 
    load_gene_annotations(
        gene.types=
            c(
                'lincRNA',
                'snRNA',
                'miRNA',
                'snoRNA',
                # 'transcribed_unprocessed_pseudogene',
                # 'processed_transcript',
                # 'transcribed_processed_pseudogene',
                'protein_coding'
            )
    ) %>% 
    standardize_data_cols(skip.resolution=TRUE)
# Now combine all annotations together?
all.functional.loci.gene.associations <- 
    check_cached_results(
        results_file=ALL_CLEAN_GENE_LOCUS_ASSOCIATIONS_FILE,
        force_redo=parsed.args$force.redo,
        results_fnc=map_all_indirect_gene_locus_associations_with_gene_metadata,
        gene.positions.df=gene.positions.df,
        gene.locus.association.dfs=
            list(
                enhancers.df,
                # eQTLs.df,
                # TFs.df,
                TFBS.df
            )
    )
# all.functional.loci.gene.associations
# all.functional.loci.gene.associations %>% count(association.type, association.subtype, association.source)

# all.functional.loci.gene.associations %>% head(2) %>% t()
