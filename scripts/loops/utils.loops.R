# library('FitHiC')
library(stringi)
library(furrr)
library(idr2d)

###################################################
# cooltools dots
###################################################
load_cooltools_dots <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE
    )
}

list_all_cooltools_dots_results <- function(){
    LOOPS_DIR %>% 
    parse_results_filelist(suffix='-dots.tsv') %>%
    filter(method == 'cooltools') %>% 
    get_info_from_MatrixIDs(keep_id=FALSE)
}

load_all_cooltools_dots <- function(){
    list_all_cooltools_dots_results() %>% 
    mutate(
        loops=
            # pmap(
            future_pmap(
                .,
                load_cooltools_dots,
                .progress=TRUE
            )
    ) %>%
        # {.} -> tmp; tmp
        # tmp %>% select(SampleID, loops)
    unnest(loops) %>% 
    dplyr::rename(chr=chrom1) %>% 
    rename_with(
        ~ stri_replace_all(
            .x, 
            regex='la_exp.([a-z]+).(value|qval)',
            '$2.$1'
        )
    ) %>% 
    select(
        -c(
            method,
            cstart1,
            cstart2,
            region,
            region1,
            region2,
            end1,
            end2,
            chrom2,
            filepath
        )
    )
}

post_process_cooltools_dots_results <- function(results.df) {
    results.df %>%
    dplyr::rename(
        'anchor.left'=start1,
        'anchor.right'=start2
    ) %>% 
    pivot_longer(
        c(
            starts_with('value.'),
            starts_with('qval.')
        ),
        names_to='statistic',
        values_to='value'
    ) %>% 
    separate_wider_delim(
        statistic,
        delim='.',
        names=c('statistic', 'kernel')
    ) %>% 
    pivot_wider(
        names_from=statistic,
        values_from=value
    ) %>%
    dplyr::rename('enrichment'=value) %>% 
    mutate(
        SampleID=str_replace_all(SampleID, '.Merged.Merged', ''),
        log10.qval=-log10(qval),
        length=anchor.right - anchor.left
    ) %>% 
    dplyr::select(
        c(
            "type",
            "weight",
            "resolution",
            # "Edit",
            # "Celltype",
            # "Genotype",
            # "CloneID",
            # "TechRepID",
            # "ReadFilter",
            # "isMerged",
            "SampleID",
            "chr",
            "anchor.left",
            "anchor.right",
            "count",
            # "c_label",
            # "c_size",
            "length",
            "kernel",
            "enrichment",
            "log10.qval"
        )
    )
}

filter_loop_results <- function(
    results.df,
    q.thresh=0.1){
    results.df %>% 
    filter(kernel == 'donut') %>% 
    filter(weight == 'balanced') %>% 
    filter(log10.qval <= -log10(q.thresh))
}

###################################################
# IDR2D Analysis
###################################################
tidy_IDR2D_sided_results <- function(
    results.obj,
    metric_colname,
    ...){
    results.obj %>% 
    as_tibble() %>% 
    mutate(
        # idr=ifelse(is.na(idr), -1, idr), # to help identify rep-exclusive loops
        diff.value=value - rep_value,
        diff.rank=rank - rep_rank
    ) %>% 
    rename(
        'chr'=chr_a,
        'anchor.left'=start_a,
        'anchor.right'=start_b,
        'IDR'=idr
    ) %>% 
    select(
        chr, anchor.left, anchor.right,
        diff.value, diff.rank,
        IDR
    )
}

tidy_IDR2D_results <- function(
    results,
    metric_colname,
    ...){
    # all loops from rep1
    reproducible.loops.P1 <- 
        tidy_IDR2D_sided_results(
            results$rep1_df,
            metric_colname
        )
    # all loops from rep2
    reproducible.loops.P2 <- 
        tidy_IDR2D_sided_results(
            results$rep2_df,
            metric_colname
        ) %>%
        # consistent sign so -ve => rep1 > rep2 for both sets of results
        mutate(
            across(
                starts_with('diff.'),
                ~ -.x
            )
        )
    # combine all loop results from both replicates
    # loops detected in both
    bind_rows(
        reproducible.loops.P1 %>% filter(!is.na(IDR)),
        reproducible.loops.P2 %>% filter(!is.na(IDR)),
    ) %>% 
    # loops detected in exactly one of the replicates
    bind_rows(
        reproducible.loops.P1 %>% 
            filter(is.na(IDR)) %>%
            add_column(loop.type='P1.only'),
        reproducible.loops.P2 %>%
            filter(is.na(IDR)) %>%
            add_column(loop.type='P2.only'),
    ) %>%
    distinct(pick(-c('loop.type')), .keep_all=TRUE) %>% 
    # create column indicating which loops are reproduced between conditions
    mutate(
        loop.type=
            case_when(
                is.na(loop.type) & !is.na(IDR) ~ 'detected.in.both',
                loop.type == 'P1.only'         ~ loop.type,
                loop.type == 'P2.only'         ~ loop.type,
                TRUE                           ~ NA
            )
    )
}

run_IDR2D_analysis <- function(
    loops.P1,
    loops.P2,
    metric_colname,
    value_transformation,
    ambiguity_resolution_method,
    max_gap,
    ...){
    # paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row_index]]", collapse='; ')
    loops.P1 <- 
        loops.P1 %>% 
        mutate(across(c(chr.A, chr.B), as.character)) %>% 
        mutate(across(c(start.A, start.B, end.A, end.B), as.integer)) %>% 
        select(
            chr.A, start.A, end.A,
            chr.B, start.B, end.B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    loops.P2 <-  
        loops.P2 %>% 
        mutate(across(c(chr.A, chr.B), as.character)) %>% 
        mutate(across(c(start.A, start.B, end.A, end.B), as.integer)) %>% 
        select(
            chr.A, start.A, end.A,
            chr.B, start.B, end.B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    # run IDR2D to define replicable loops between conditions
    estimate_idr2d(
        loops.P1,
        loops.P2, 
        value_transformation=value_transformation,
        ambiguity_resolution_method=ambiguity_resolution_method,
        max_gap=ifelse(max_gap < 1, -1L, max_gap)
    ) %>% 
    # tidy up results into nice tabular format
    tidy_IDR2D_results(
        loops.P1,
        loops.P2,
        metric_colname
    )
}

run_all_IDR2D_analysis <- function(
    nested.loops.df,
    hyper.params.df,
    force.redo,
    sample.group.comparisons,
    pair_grouping_cols,
    SampleID.fields,
    sampleID_col='SampleID',
    suffixes=c('.P1', '.P2'),
    ...){
    # force.redo=parsed.args$force.redo; sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS %>% rename('SampleID.P1'=Sample.Group.Numerator, 'SampleID.P2'=Sample.Group.Denominator); pair_grouping_cols=c('weight', 'resolution', 'kernel', 'chr'); SampleID.fields=c(NA, 'Celltype', 'Genotype'); sampleID_col='SampleID'; suffixes=c('.P1', '.P2')
    # list + format metadata for all specified sample groups to compare
    nested.loops.df %>% 
    enumerate_pairwise_comparisons(
        sample.group.comparisons=sample.group.comparisons,
        pair_grouping_cols=pair_grouping_cols,
        sampleID_col=sampleID_col,
        suffixes=suffixes,
        SampleID.fields=SampleID.fields,
        include_merged_col=FALSE
    ) %>% 
    # Evaluate all comparisons for all combinations of specified parameters
    cross_join(hyper.params.df) %>% 
    mutate(
        max_gap=resolution * max_gap_bins,
        output_dir=
            file.path(
                LOOPS_IDR2D_DIR,
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('kernel_{kernel}'),
                glue('metric_{metric_colname}'),
                glue('resolve.method_{ambiguity_resolution_method}'),
                glue('max.gap_{max_gap}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.P1}_vs_{SampleID.P2}-IDR2D.tsv')
            )
    ) %>% 
    arrange(desc(chr)) %>% 
    {
        if (!force.redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            mutate(
                ., 
                idr.params=glue('{metric_colname}-{value_transformation}-{ambiguity_resolution_method}'),
                comparison=glue('{SampleID.P1}-{SampleID.P2}')
            ) %>% 
            dplyr::count(
                resolution, max_gap,
                idr.params,
                comparison
                # .,
                # metric_colname, value_transformation, ambiguity_resolution_method,
                # resolution, max_gap,
                # SampleID.P1, SampleID.P2
            )
        )
    } %>%
        # {.} -> tmp; tmp
        # tmp %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        force_redo=force.redo,
        return_data=FALSE,
        results_fnc=run_IDR2D_analysis,
        .progress=TRUE
    )
}

load_IDR2D_results <- function(filepath, ...){
    filepath %>%
    read_tsv(show_col_types=FALSE)
}

list_all_IDR2D_results <- function(
    suffixes=c('.P1', '.P2'),
    ...){
    LOOPS_IDR2D_DIR %>% 
    parse_results_filelist(
        suffix='-IDR2D.tsv',
        filename.column.name='Sample.Group.Pair'
    ) %>% 
    separate_wider_delim(
        Sample.Group.Pair,
        delim='_vs_',
        names=paste0('SampleID', suffixes)
    )
}

load_all_IDR2D_results <- function(){
    list_all_IDR2D_results() %>% 
    mutate(
        idr2d=
            # pmap(
            future_pmap(
                .,
                load_IDR2D_results,
                .progress=TRUE
            )
    ) %>%
    unnest(idr2d) %>% 
    select(-c(filepath))
}

post_process_IDR2D_results <- function(results.df){
    results.df %>% 
    mutate(
        comparison=glue('{SampleID.P1} vs {SampleID.P2}'),
        max.gap.bins.int=max.gap / resolution,
        max.gap=
            fct_reorder(
                glue('max gap = {scale_numbers(max.gap, force_chr=TRUE)}'),
                max.gap.bins.int
            ),
        max.gap.bins=
            fct_reorder(
                glue('max gap = {max.gap.bins.int} bins'),
                max.gap.bins.int
            ),
        is.loop.shared=
            case_when(
                loop.type == 'P1.only' ~ 'P1.only',
                loop.type == 'P2.only' ~ 'P2.only',
                IDR <= 0.1             ~ 'IDR < 0.1',
                IDR <= 1               ~ 'Irreproducible',
                TRUE                   ~ NA
            ) %>%
            factor(
                levels=
                    c(
                        'IDR < 0.1',
                        'Irreproducible',
                        'P1.only',
                        'P2.only'
                    )
            )
    ) %>%
    select(
        -c(
            loop.type,
            region,
            weight,
            kernel,
            metric, 
            resolve.method
        )
    )
}

filter_loop_IDR2D_results <- function(results.df){
    results.df %>% 
    filter(weight == 'balanced') %>% 
    filter(kernel == 'donut') %>% 
    filter(metric == 'log10.qval') %>% 
    filter(resolve.method == 'value')
}

###################################################
# Expression Analysis
###################################################
join_expr_and_IDR2D_results <- function(
    idr2d.results.df,
    expression.df,
    ...){
    # samples.avail <- expression.df %>% colnames() %>% grep('16p.', ., value=TRUE)
    samples.avail <- unique(expression.df$SampleID)
    idr2d.results.df %>% 
    filter(
           SampleID.P1 %in% samples.avail,
           SampleID.P2 %in% samples.avail
    ) %>% 
    inner_join(
        expression.df,
        by=
            join_by(
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    )
}

calc_expr_loop_ztest <- function(
    idr2d.results.df,
    expression.df,
    ...){
    idr2d.results.df %>% 
    inner_join(
        expression.df,
        by=
            join_by(
                SampleID.P1 == SampleID,
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    ) %>% 
    inner_join(
        expression.df,
        suffix=c('.P1', '.P2'),
        by=
            join_by(
                SampleID.P2 == SampleID,
                chr,
                start, end,
                symbol, EnsemblID
            )
    ) %>%
    # for each gene compute pvalue if mean expression is different between conditions
    add_column(
        n.rna.replicates.P1=6,
        n.rna.replicates.P2=6
    ) %>% 
    mutate(
        TPM.se.P1=TPM.sd.P1**2 / n.rna.replicates.P1,
        TPM.se.P2=TPM.sd.P2**2 / n.rna.replicates.P2,
        expr.Z=(TPM.mean.P2 - TPM.mean.P1) / sqrt(TPM.se.P1 + TPM.se.P2),
        expr.p=2 * (1 - pnorm(abs(expr.Z)))
    ) %>%
    # adjust p-values genome-wide
    group_by(
        weight, resolution, kernel,
        resolve.method, metric,
        SampleID.P1, SampleID.P2
    ) %>% 
    mutate(expr.p.adjust=p.adjust(expr.p, method='BH')) %>% 
    ungroup() %>% 
    select(
        -c(
            n.rna.replicates.P1, n.rna.replicates.P2,
            TPM.se.P1, TPM.se.P2,
            # TPM.sd.P1, TPM.sd.P2,
            # TPM.mean.P1, TPM.mean.P2,
            expr.Z, expr.p
        )
    )
}

count_ccres_per_loop <- function(
    ccres.df=NULL,
    idr2d.df=NULL, 
    regions.df=NULL){
    # Specific sub-regions of each chrosome to quantify stats for separately
    if (is.null(regions.df)) {
        regions.df <- 
            GENOMIC_REGIONS %>% 
            filter(region != 'chr16') %>% 
            dplyr::rename_with(~ str_remove(.x, 'region.')) %>%
            select(region, chr, start, end)
    }
    # cis-regulatory element annotations (cCREs)
    if (is.null(ccres.df)) {
        # All cCREs are 150bp - 350bp long
        # ccres.df %>% mutate(dist=end - start) %>% group_by(cCRE.type) %>% summarize(as_tibble_row(summary(dist)))
        ccres.df <- 
            load_encode_ccres() %>%
            left_join(
                regions.df,
                suffix=c('', '.y'),
                by=
                    join_by(
                        chr,
                        within(
                            x$start, x$end,
                            y$start, y$end
                        )
                    )
            ) %>% 
            mutate(region=ifelse(is.na(region), chr, region)) %>% 
            select(-c(start.y, end.y))
    }
    # loop reproducibility results
    if (is.null(idr2d.df)) {
        idr2d.df <- 
            check_cached_results(
                results_file=LOOPS_IDR2D_RESULTS_FILE,
                # force_redo=TRUE,
                results_fnc=load_all_IDR2D_results
            ) %>% 
            post_process_IDR2D_results() %>% 
            standardize_data_cols() %>% 
            left_join(
                regions.df,
                by=
                    join_by(
                        chr,
                        within(
                            x$anchor.left, x$anchor.right,
                            y$start, y$end
                        )
                    )
            ) %>% 
            mutate(region=ifelse(is.na(region), chr, region)) %>% 
            select(
                -c(
                    max.gap,
                    SampleID.P1,
                    SampleID.P2,
                    start,
                    end
                )
            )
    }
    # total number of cCREs per category per region
    ccre.totals.df <- 
        ccres.df %>% 
        count(
            chr, region,
            cCRE.type,
            name='total.cCREs'
        )
    # total number of loops per category per region
    loop.totals.df <- 
        idr2d.df %>% 
        count(
            resolution, max.gap.bins,
            comparison,
            chr, region,
            is.loop.shared,
            name='total.loops'
        )
    # loop.totals.df
    # ccre.totals.df
    idr2d.df %>%
        filter(max.gap.bins.int == 5) %>% 
    # nest so one set of loop annotations per row
    nest(
        loops=
            c(
                # is.loop.shared,
                diff.value,
                diff.rank,
                IDR,
                anchor.left,
                anchor.right
            )
    ) %>% 
    # match loop sets and cCRE sets by region
    left_join(
        # nest so one set of cCREs annotations per row
        ccres.df %>%
        nest(
            cCREs=
                c(
                    # cCRE.type,
                    start,
                    end,
                    cCREID
                )
        ),
        by=
            join_by(
                chr,
                region
            )
    ) %>% 
    # For each region, compute all overlaps between any loop and any cCRE
    rowwise() %>% 
    mutate(
        overlaps=
            inner_join(
                loops,
                cCREs,
                by=
                    join_by(
                        within(
                            y$start, y$end,
                            x$anchor.left, x$anchor.right
                        )
                    )
            ) %>%
            list()
    ) %>% 
    ungroup() %>% 
        {.} -> tmp; tmp
        tmp %>% 
            rowwise() %>% mutate(across(c(loops, cCREs, overlaps), ~ nrow(.x))) %>% ungroup() %>%
            dplyr::rename('max.gap'=max.gap.bins.int) %>% 
            select(
                   resolution, max.gap,
                   comparison,
                   region,
                   is.loop.shared, loops,
                   cCRE.type, cCREs,
                   overlaps
            )
        tmp %>% 
            rowwise() %>% mutate(across(c(loops, cCREs, overlaps), ~ nrow(.x))) %>% 
            ungroup() %>% summarize(across(c(loops, cCREs, overlaps), sum))
}

