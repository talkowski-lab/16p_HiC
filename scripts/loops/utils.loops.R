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

post_process_cooltools_dots_results <- function(results.df){
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
            # "type",
            "weight",
            "resolution",
            "Edit",
            "Celltype",
            "Genotype",
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

###################################################
# Analysis
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
    distinct() %>% 
    # loops detected in exactly one of the replicates
    bind_rows(
        reproducible.loops.P1 %>% 
            filter(is.na(IDR)) %>%
            add_column(loop.type='P1.only'),
        reproducible.loops.P2 %>%
            filter(is.na(IDR)) %>%
            add_column(loop.type='P2.only'),
    ) %>%
    # create column indicating which loops are reproduced between conditions
    mutate(
        loop.type=
            case_when(
                is.na(loop.type) & IDR < 0.01  ~ 'IDR < 0.01',
                is.na(loop.type) & IDR < 0.05  ~ 'IDR < 0.05',
                is.na(loop.type) & IDR < 0.1   ~ 'IDR < 0.1 ',
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
    ...){
    # paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[1]]", collapse=';')
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
        # ...,
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
    # all.loops.df=loops.df; sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS %>% rename('SampleID.P1'=Sample.Group.Numerator, 'SampleID.P2'=Sample.Group.Denominator); pair_grouping_cols=c('weight', 'resolution', 'kernel', 'chr'); SampleID.fields=c(NA, 'Celltype', 'Genotype'); sampleID_col='SampleID'; suffixes=c('.P1', '.P2')
    # list + format metadata for all specified sample groups to compare
    nested.loops.df %>% 
    enumerate_pairwise_comparisons(
        sample.group.comparisons=sample.group.comparisons,
        pair_grouping_cols=pair_grouping_cols,
        sampleID_col=sampleID_col,
        suffixes=suffixes,
        SampleID.fields=SampleID.fields,
        include_merged_col=FALSE,
        keep_id=FALSE
    ) %>% 
    # Evaluate all comparisons for all combinations of specified parameters
    cross_join(hyper.params.df) %>% 
    mutate(
        output_dir=
            file.path(
                LOOPS_IDR2D_DIR,
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('kernel_{kernel}'),
                glue('metric_{metric_colname}'),
                glue('resolve.method_{ambiguity_resolution_method}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.P1}_vs_{SampleID.P2}-IDR2D.tsv')
            )
    ) %>% 
    arrange(desc(chr)) %>% 
        # {.} -> tmp; tmp
        # tmp %>% head(3) %>% 
            # pmap(
    future_pmap(
        .l=.,
        .f=check_cached_results,
        force_redo=force.redo,
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
    mutate(
        loop.type=
            factor(
                loop.type, 
                levels=
                    c(
                        'IDR < 0.01',
                        'IDR < 0.05',
                        'IDR < 0.1',
                        'detected.in.both',
                        'P1.only',
                        'P2.only'
                    )
            )
    ) %>% 
    select(
        -c(
            region,
            filepath
        )
    )
}

###################################################
# Plotting
###################################################
plot_upset <- function(
    plot.df,
    category_prefix,
    title.str,
    ...){
    upset(
        plot.df,
        plot.df %>%
            dplyr::select(starts_with(category_prefix)) %>%
            colnames(),
        width_ratio=0.3,
        mode='exclusive_intersection',
        name=category_prefix,
        # name=category_col,
        # labeller=function(x) str_remove(x, category_prefix),
        annotations=
            list(
                'Chrs'=
                    (
                        ggplot(mapping=aes(fill=chr))
                        + geom_bar(stat='count', position='fill')
                        + scale_y_continuous(labels=scales::percent_format())
                        + ylab('Chrs')
                    )
            ),
        set_sizes=
            (
                upset_set_size(
                    position='right',
                    geom=
                        geom_bar(
                            aes(fill=chr, x=group),
                            width=0.8
                        )
                ) +
                make_ggtheme(axis.text.x=element_text(angle=45, hjust=1))
            ),
        guides='over' # moves legends over the set sizes
    ) +
    ggtitle(title.str)
}

plot_pairs <- function(
    freq.df,
    cols_pattern,
    text.size=7,
    ...){
    cols_to_plot <-
        freq.df %>%
        colnames() %>%
        str_detect(cols_pattern) %>%
        which()
    ggpairs(
        freq.df,
        columns=cols_to_plot,
        ...
    ) +
    theme(
        axis.text.x=element_text(size=text.size, angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=text.size),
        strip.text=element_text(size=text.size)
    )
}



