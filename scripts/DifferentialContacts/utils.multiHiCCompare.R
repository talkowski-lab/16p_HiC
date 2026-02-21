###################################################
# Depdendencies
###################################################
# library(tidyverse)
# library(magrittr)
library(glue)
library(multiHiCcompare)
library(BiocParallel)
# library(ggplot2)
library(viridis)
library(cowplot)
library(gtable)
library(purrr)
library(furrr)
library(ComplexUpset)
# library(hictkR)

###################################################
# Generate resutls
###################################################
handle_covariates <- function(
    samples.df,
    covariates.df,
    effect.col='Sample.Group',
    sampleID.col='SampleID'){
    # Immediate check 
    if (is.null(covariates.df)) {
        covariates <- NULL
        design.matrix <- NULL
        message('No covaraites to use')
    } else {
        # List all covariates for this set of samples
        covariates <- 
            samples.df %>% 
            select(all_of(c(sampleID.col, effect.col))) %>% 
            left_join(
                covariates.df,
                by=sampleID.col
            ) %>%
            select(-all_of(sampleID.col))
        # Check that batch variables are not uniform 
        covariate_level_sizes <- 
            covariates %>% 
            select(-c(all_of(effect.col))) %>% 
            pivot_longer(everything(), names_to='covariate', values_to='value') %>%
            distinct() %>% 
            count(covariate) %>% 
            deframe()
        # If all covariates are uniform, return null, no design matrix
        if (all(covariate_level_sizes == 1)) {
            covariates <- NULL
            design.matrix <- NULL
            message('Covariates are uniform, ignoring')
        } else {
            # List all uniform covariates (uninformative)
            non.uniform.covariates <- 
                covariate_level_sizes[covariate_level_sizes != 1] %>% names()
            # Make formula for GLM
            covariate.names <- 
                covariates %>% 
                select(-c(all_of(effect.col))) %>% 
                colnames() %>%
                {.[. %in% non.uniform.covariates]}
            contrast <- 
                c(effect.col, covariate.names) %>%
                paste(collapse=" + ") %>% 
                sprintf('~ %s', .) %>% 
                formula()
            design.matrix <- model.matrix(contrast, covariates)
        }
    }
    # return final design matrix
    return(
        list(
            covariates=covariates,
            design.matrix=design.matrix
        )
    )
}

run_multiHiCCompare <- function(
    samples.df,
    sample_group_priority_fnc,
    zero.p,
    A.min,
    resolution,
    range1,
    range2,
    remove.regions,
    md_plot_file,
    covariates.df=NULL,
    frac.cutoff=0.8,
    effect.col='Sample.Group',
    p.method='fdr',
    ...){
    # row_index=2; sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC; samples.df=tmp$samples.df[[row_index]]; zero.p=tmp$zero.p[[row_index]]; A.min=tmp$A.min[[row_index]]; resolution=tmp$resolution[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]; md_plot_file=tmp$md_plot_file[[row_index]]; remove.regions=hg38_cyto; covariates.df=NULL; frac.cutoff=0.8; effect.col='Sample.Group'; p.method='fdr';
    # Handle covariates if specified
    design.info <- 
        handle_covariates(
            samples.df,
            covariates.df,
            effect.col
        )
    # Get sample groups, ensure consistent num/denom for fc estimates
    # function with more priority (smaller number) will be numerator
    # should be no ties, but ties are broken alphabetically
    samples.df <- 
        samples.df %>%
        mutate(
            group.priority=sample_group_priority_fnc(!!sym(effect.col)),
            !!effect.col :=
                fct_reorder(
                    !!sym(effect.col),
                    group.priority,
                    .desc=TRUE
                )
        ) %>%
        arrange(!!effect.col)
    # Load all contacts for samples + regions
    samples.contacts <- 
        samples.df %>%
        # select(-c(resolution)) %>% 
        pmap(
            .l=.,
            .f=load_mcool_file,
            resolution=resolution,
            range1=range1,
            range2=range2,
            normalization="NONE",
            .progress=FALSE
        )
    message('loaded contact matrices')
    # Get all bin pairs that are detected in > frac.cutoff fraction of samples e.g. 
    # frac.cutoff=0.8 -> only test contacts detected in > 80% of all samples
    common.bin.pairs <- 
        samples.contacts %>%
        bind_rows(.id='index') %>%
        select(chr, range1, range2) %>% 
        count(chr, range1, range2) %>% 
        filter(n >= max(n) * frac.cutoff) %>%
        select(chr, range1, range2)
    message(glue('Testing {nrow(common.bin.pairs)} bin-pairs for DAC'))
    # Now only subset to commonly found contacts 
    samples.contacts <- 
        samples.contacts %>%
        lapply(
            function(df) {
                inner_join(
                    df,
                    common.bin.pairs,
                    by=join_by(chr, range1, range2)
                )
            }
        )
    message(glue('subset contacts to those detected in {frac.cutoff * 100}% of contacts'))
    # Make experiment object with relevant data+parameters
    make_hicexp(
        data_list=samples.contacts,
        groups=samples.df %>% pull(!!effect.col),
        covariates=design.info$covariates,
        zero.p=zero.p,
        A.min=A.min,
        remove.regions=remove.regions,
        remove_zeroes=FALSE,
        filter=TRUE
    ) %T>% 
    {message('Made experiment object')} %>% 
    # Normalize hic data, use cyclic loess with automatically calculated span
    fastlo(
        verbose=TRUE,
        parallel=TRUE,
        span=NA
    ) %T>% 
    # Plot normalized IFs for all sample pairs
    {
        message('Normalized contacts with fast cyclic loess')
        dir.create(
            dirname(md_plot_file),
            showWarnings=FALSE,
            recursive=TRUE
        )
        pdf(
            md_plot_file,
            height=nrow(.@metadata) * 2,
            width=6
        )
        MD_hicexp(
            .,
            pcol=2
        )
        dev.off()
    } %>%
    # Preform differential testing on contacts (bin-pairs)
    # Handle covariate information
    {
        if (is.null(design.info$design.matrix)) {
            hic_exactTest(
                .,
                p.method=p.method,
                parallel=TRUE
            )
        } else {
            hic_glm(
                .,
                design=design.info$design.matrix,
                coef=2,  # Sample.Group (only correct if 2 levels), 1 would be intercept
                method="QLFTest",
                p.method=p.method,
                parallel=TRUE
            )
        }
    } %>% 
    # Get differential results
    results() %>% 
    as_tibble() %>%
    arrange(p.adj)
}

run_all_multiHiCCompare <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    # remove.regions,
    force_redo=FALSE,
    covariates.df=NULL,
    chromosomes=CHROMOSOMES,
    group1_colname='Sample.Group.P1',
    group2_colname='Sample.Group.P2',
    ...){
    # sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC; force_redo=FALSE; covariates.df=NULL; chromosomes=c('chr15', 'chr16'); group1_colname='Sample.Group.P1'; group2_colname='Sample.Group.P2'
    comparisons.df %>% 
    # for each Sample.Group list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    # creates the columns Sample.Group.Numerator, Sample.Group.Denominator[
    set_foldchange_direction_as_factor(
        sample_group_priority_fnc=sample_group_priority_fnc,
        group1_colname=group1_colname,
        group2_colname=group2_colname,
    ) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                MULTIHICCOMPARE_DIR,
                'results',
                # glue('merged_{isMerged}'),
                glue('zero.p_{zero.p}'),
                glue('A.min_{A.min}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                # glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        md_plot_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-MD.plot.pdf')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-multiHiCCompare.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !(file.exists(results_file)))
        } else{
            .
        }
    } %>% 
    arrange(desc(chr)) %>% 
    arrange(desc(resolution)) %T>% 
    {
        message('Generating the following results files')
        print(
            dplyr::count(
                .,
                zero.p,
                A.min,
                resolution,
                Sample.Group.Numerator,
                Sample.Group.Denominator
            ) %>%
            rename_with(~ str_remove(., 'Sample.Group.'))
        )
    } %>%
        # {.} -> tmp
        # tmp %>% head(2) %>% 
    future_pmap(
    # pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_multiHiCCompare,
                    sample_group_priority_fnc=sample_group_priority_fnc,
                    covariates.df=covariates.df,
                    # all columns also passed as input to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load resutls
###################################################
load_and_correct_multiHiCCompare_results <- function(
    filepaths,
    nom.threshold,
    fdr.threshold,
    gw.fdr.threshold,
    ...){
    filepaths %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    bind_rows() %>%
    mutate(p.adj.gw=p.adjust(p.value, method='BH')) %>% 
    filter(
        p.adj.gw < gw.fdr.threshold,
        p.adj    < fdr.threshold,
        p.value  < nom.threshold
    )
}

list_all_multiHiCCompare_results <- function(
    resolutions=NULL,
    sample.group.comparisons=NULL,
    file_suffix='-multiHiCCompare.tsv',
    ...){
    # resolutions=NULL; sample.group.comparisons=NULL; file_suffix='-multiHiCCompare.tsv'
    if (is.null(sample.group.comparisons)){
        sample.group.comparison.colnames <- c('Sample.Group.Numerator','Sample.Group.Denominator')
    } else {
        sample.group.comparison.colnames <- colnames(sample.group.comparisons)
    }
    # Load all results
    parse_results_filelist(
        input_dir=file.path(MULTIHICCOMPARE_DIR, 'results'),
        suffix=file_suffix,
        filename.column.name='pair.name',
        param_delim='_',
    ) %>%
    # Split title into pair of groups ordered by numerator/denominator
    mutate(pair.name=str_remove(pair.name, file_suffix)) %>% 
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=sample.group.comparison.colnames
    ) %>% 
    mutate(
        Sample.Group=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}'),
        resolution=scale_numbers(resolution, force_numeric=TRUE)
    ) %>% 
    # Filter relevant results sets
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(sample.group.comparisons)) {
            inner_join(
                .,
                sample.group.comparisons,
                by=sample.group.comparison.colnames
            )
        } else {
            .
        }
    }
}

load_all_multiHiCCompare_results <- function(
    resolutions=NULL,
    sample.group.comparisons=NULL,
    gw.fdr.threshold=1,
    fdr.threshold=1,
    nom.threshold=1,
    ...){
    # sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC; sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS; gw.fdr.threshold=0.1; fdr.threshold=0.1; nom.threshold=0.05; resolutions=NULL
    list_all_multiHiCCompare_results(
        resolutions=resolutions,
        sample.group.comparisons=sample.group.comparisons
    ) %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    group_by(across(-c(filepath, region))) %>% 
    summarize(filepaths=list(c(filepath))) %>% 
    ungroup() %>% 
    mutate(
        results=
            # future_pmap(
            pmap(
                .l=.,
                # correct adjusted pvalues genome-wide
                .f=load_and_correct_multiHiCCompare_results,
                # Only load results with specific params
                gw.fdr.threshold=gw.fdr.threshold,
                fdr.threshold=fdr.threshold,
                nom.threshold=nom.threshold,
                .progress=TRUE
            )
    ) %>% 
    unnest(results) %>% 
    select(-c(D, filepaths)) %>%
    mutate(
    # calculate log of all pvalues + add columns
        across(
            starts_with('p.'),
            .f=~ -log10(.x),
            .names='log.{.col}'
        )
    ) %>% 
    mutate(
        chr=
            chr %>% 
            rename_chrs(to_label=TRUE) %>% 
            factor(levels=CHROMOSOMES)
    ) %>% 
    mutate(
        tidy.metadata=
            tidy_pair_metadata(
                sampleID.pairs.df=
                    select(
                        .data=., 
                        all_of(c('Sample.Group.Numerator', 'Sample.Group.Denominator'))
                    ),
                suffixes=c('.Numerator', '.Denominator'),
                SampleID.fields=c('Edit', 'Celltype', 'Genotype'),
                # ...
            )
    ) %>%
    unnest(tidy.metadata)
}

post_process_multiHiCCompare_results <- function(
    results.df,
    chromosomes=CHROMOSOMES){
    # CHROMOSOMES
    # CHROMOSOMES[c(1,5,10,16,22)]
    # CHROMOSOMES[c(10,16,22)]
    results.df %>% 
    # filter results 
    filter(chr %in% chromosomes) %>% 
    filter(A.min == 5) %>% 
    filter(zero.p == 0.8) %>% 
    select(-c(A.min, zero.p)) %>% 
    mutate(distance.bp=region2 - region1) %>% 
    unite(
        'bin.pair.idx',
        chr, region1, region2,
        sep='#',
        remove=FALSE
    )
}

load_correct_count_multiHiCCompare_results <- function(
    filepaths,
    ...){
    filepaths %>% 
    load_and_correct_multiHiCCompare_results(
        nom.threshold=1,
        fdr.threshold=1,
        gw.fdr.threshold=1
    ) %>% 
    # calculate log of all pvalues + add columns
    mutate(
        across(
            starts_with('p.'),
            .f=~ -log10(.x),
            .names='log.{.col}'
        )
    ) %>% 
   mutate(
        'sig.lvl.p.adj.gw < 1e-15'=p.adj.gw <  1e-15,
        'sig.lvl.p.adj.gw < 1e-10'=p.adj.gw <  1e-10,
        'sig.lvl.p.adj.gw < 1e-5' =p.adj.gw <  1e-5,
        'sig.lvl.p.adj.gw < 0.001'=p.adj.gw <  1e-3,
        'sig.lvl.p.adj.gw < 0.1'  =p.adj.gw <  0.1,
        'sig.lvl.N.S.'            =p.adj.gw >= 0.1
    ) %>% 
    pivot_longer(
        starts_with('sig.lvl.'),
        names_to='sig.lvl',
        names_prefix='sig.lvl.',
        values_to='meet.sig.lvl'
    ) %>% 
    filter(meet.sig.lvl) %>% 
    count(
        # A.min, zero.p, merged,
        # resolution, 
        # Sample.Group, 
        # Celltype,
        chr,
        sig.lvl,
        name='nDACs'
    ) %>%
   mutate(
        sig.lvl=
            factor(
                sig.lvl,
                levels=
                    c(
                        'N.S.',
                        'p.adj.gw < 0.1',
                        'p.adj.gw < 0.001',
                        'p.adj.gw < 1e-5',
                        'p.adj.gw < 1e-10',
                        'p.adj.gw < 1e-15'
                    )
            )
    )
}

count_contacts_by_significance <- function(
    resolutions=NULL,
    sample.group.comparisons=NULL,
    ...){
    list_all_multiHiCCompare_results(
        resolutions=resolutions,
        sample.group.comparisons=sample.group.comparisons
    ) %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    group_by(across(-c(filepath, region))) %>% 
    summarize(filepaths=list(c(filepath))) %>% 
    ungroup() %>% 
    mutate(
        results=
            # future_pmap(
            pmap(
                .l=.,
                # correct adjusted pvalues genome-wide
                .f=load_correct_count_multiHiCCompare_results,
                .progress=TRUE
            )
    ) %>% 
    unnest(results) %>% 
    select(-c(filepaths)) %>%
    mutate(
        chr=
            chr %>% 
            rename_chrs(to_label=TRUE) %>% 
            factor(levels=CHROMOSOMES)
    ) %>% 
    mutate(
        tidy.metadata=
            tidy_pair_metadata(
                sampleID.pairs.df=
                    select(
                        .data=., 
                        all_of(c('Sample.Group.Numerator', 'Sample.Group.Denominator'))
                    ),
                suffixes=c('.Numerator', '.Denominator'),
                SampleID.fields=c('Edit', 'Celltype', 'Genotype'),
                # ...
            )
    ) %>%
    unnest(tidy.metadata)
}

###################################################
# Plotting
###################################################
plot_upset <- function(
    plot.df,
    make.binary=FALSE,
    category_col='Sample.Group',
    title.str='Common DACs across Comparisons',
    ...){
    category_prefix <- fixed(glue('{category_col}.'))
    if (make.binary) {
        plot.df <-
            plot.df %>%
            add_column(is.category=TRUE) %>%
            pivot_wider(
                names_from=category_col,
                names_prefix=category_prefix,
                values_from=is.category,
                values_fill=FALSE
            )
    } 
    upset(
        plot.df,
        plot.df %>%
            dplyr::select(starts_with(category_prefix)) %>%
            colnames(),
        width_ratio=0.3,
        mode='exclusive_intersection',
        name=category_col,
        labeller=function(x) str_remove(x, category_prefix),
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

n_plot_fnc <- 
    partial(
        make_nested_plot_tabs,
        max.header.lvl=3,
        plot.fnc=plot_heatmap,
        # x.var='sig.lvl',
        fill.var='nDACs', 
        label.var='nDACs',
        # legend.position='top',
        # axis.text.x=element_text(angle=45, hjust=1),
        # legend.text=element_text(angle=35, hjust=1),
        fill.scale.mode='m',
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
    )

dist_plot_fnc <- 
    partial(
        make_nested_plot_tabs,
        plot.fnc=plot_boxplot,
        # group.cols=c('Celltype', 'resolution'),
        max.header.lvl=3,
        # x.var='Sample.Group',
        # x.scale.mode='discrete',
        y.var='value',
        facet.row='statistic',
        # fill.var='Sample.Group',
        # facet.row='statistic',
        outlier.size=0.2,
        # legend.position='right',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # axis.text.x=element_text(angle=45, hjust=1),
    )

volcano_plot_fnc <- 
    partial(
        make_nested_plot_tabs,
        plot.fnc=plot_jitter,
        # group.cols=c('Celltype', 'resolution'),
        max.header.lvl=3,
        x.var='logFC',
        y.var='log.p.adj.gw',
        x.axis.label.accuracy=0.01,
        y.axis.label.accuracy=0.1,
    )

distance_plot_fnc <- 
    partial(
        make_nested_plot_tabs,
        plot.fnc=plot_jitter,
        max.header.lvl=3,
        x.var='logFC',
        y.var='distance.bp',
        y.scale.mode='mb',
        scales='free_y'
        # axis.title.y=element_blank(),
        # axis.title.x=element_blank(),
    )

manhattan_plot_fnc <- 
    partial(
        make_nested_plot_tabs,
        plot.fnc=plot_jitter,
        max.header.lvl=3,
        x.var='region.bp',
        x.scale.mode='mb',
        y.var='value',
        scale.fill.mode='mb',
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.elements=
            list(
                # geom_vline(
                #     xintercept=c(86435256, 86521792), # WAPL region on chr10
                #     linewidth=0.3,
                #     linetype='dashed',
                #     color='#619CFF'
                # ),
                # geom_vline(
                #     xintercept=c(36876769, 37066413), # NIPBL region on chr5
                #     linewidth=0.3,
                #     linetype='dashed',
                #     color='#F8766D'
                # ),
                geom_hline(
                    yintercept=0,
                    linewidth=0.05,
                    linetype='solid',
                    color='black',
                )
            )
    )

