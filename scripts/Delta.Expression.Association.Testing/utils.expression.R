###################################################
# Load Differential Expression Reuslts
###################################################
load_all_DESeq2_results <- function(force.redo=FALSE){
    check_cached_results(
        results_file=DESEQ2_RESULTS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(suffix='-DESeq2.tsv') {
                # List all results files
                DESEQ2_DATA_DIR %>% 
                list.files(
                    pattern=suffix,
                    full.names=TRUE,
                    recursive=TRUE
                ) %>% 
                tibble(filepath=.) %>%
                mutate(info=str_remove(basename(filepath), suffix)) %>% 
                filter(!grepl('iPSC', info)) %>% 
                # Tidy pairwise metadata
                separate_wider_delim(
                    info,
                    # delim='-',
                    delim='_vs_',
                    names=c('Sample.Group.Numerator', 'Sample.Group.Denominator'),
                ) %>%
                # load DEG results
                mutate(
                    results=
                        pmap(
                            list(filepath),
                            read_tsv,
                            show_col_types=FALSE
                        )
                ) %>%
                unnest(results) %>%
                # rename_with(.f=~ str_replace(.x, '^', 'DESeq2.')) %>% 
                # # run TRADEtools to define transcriptome-wide effects
                # clean up columns
                dplyr::rename('EnsemblID'=ensemblid) %>% 
                mutate(gene.length=end - start) %>% 
                select(
                    -c(
                        filepath,
                        Sample.Group.Numerator, Sample.Group.Denominator,
                        row_index,
                        # strand,
                        geneid,
                        stat
                    )
                )
            }
    ) %>% 
    standardize_data_cols()
}

###################################################
# Generate TRADE results
###################################################
run_TRADE_analysis <- function(
    deg.results,
    ...){
    print(dim(deg.results))
    # run TRADE 
    trade.results <- 
        TRADE(
            mode='univariate',
            results1=deg.results
        )
    # Significance testing results
    list(
        trade.results$significant_genes_FDR,
        trade.results$significant_genes_Bonferroni
    ) %>% 
    # Get stats + list of significant genes
    lapply(
        X=.,
        FUN=
            function(info.list) {
                # get list of significant genes
                sig.genes <- 
                    info.list %>% 
                    {.[str_detect(names(.), 'significant_gene_results_*')]} %>%
                    {.[[1]]} %>% 
                    # pull(ensemblid)
                    pull(EnsemblID)
                # get summary stats and tidy 
                info.list %>% 
                {.[!str_detect(names(.), 'significant_gene_results_*')]} %>%
                as_tibble_row() %>% 
                pivot_longer(
                    everything(),
                    names_to='stat',
                    values_to='value'
                ) %>%
                separate_wider_delim(
                    stat,
                    delim='_',
                    names=c('stat', 'is.sig', 'adjust.method')
                ) %>%
                mutate(stat=glue('{stat}.{is.sig}')) %>% 
                select(-c(is.sig)) %>% 
                pivot_wider(
                    names_from=stat,
                    names_sep='_',
                    values_from=value
                ) %>% 
                add_column('sig.genes'=list(sig.genes))
            }
    ) %>% 
    bind_rows() %>% 
    # add other stats as column
    bind_cols(
        list(
            'TWI'=trade.results$distribution_summary$transcriptome_wide_impact,
            'Effective.n.DEGs'=trade.results$distribution_summary$Me,
            'log.likelyhood'=trade.results$fit$loglik
        ) %>%
        as_tibble_row()
    ) %>%
    rename_with(~ str_replace(.x, '^', 'TRADE.'))
}

load_all_TRADE_analysis <- function(force.redo=FALSE){
    check_cached_results(
        results_file=TRADE_RESULTS_FILE,
        force_redo=force.redo,
        results_fnc=
            function() {
                load_all_DESeq2_results() %>% 
                nest(deg.results=-c(comparison, Celltype, Genotype)) %>% 
                mutate(
                    trade.df=
                        pmap(
                            .l=list(deg.results),
                             run_TRADE_analysis,
                            .progress=TRUE
                        )
                ) %>%
                unnest(trade.df) %>%
                select(-c(deg.results))
            }
    )
}

