###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(multiHiCcompare)
library(furrr)
library(BiocParallel)
library(ggplot2)
library(viridis)
library(cowplot)
library(gtable)
library(furrr)

###################################################
# Generate resutls
###################################################
set_up_sample_comparisons <- function(...){
    # get info + filepaths for all contact matrices
    load_mcool_files(
        return_metadata_only=TRUE,
        keep_metadata_columns=FALSE
    ) %>%
    # Define minimum viable resolution for each matrix
    get_min_resolution_per_matrix(as_int=TRUE, filter_res=TRUE) %>% 
    # Create groups of samples with all the same celltype+genotype, even across edits
    bind_rows(
        filter(., !isMerged) %>% 
        group_split(Genotype, Celltype, .keep=TRUE) %>%
        sapply(
            function(df) {
                ifelse(
                    length(unique(df$Edit)) > 1,
                    df,
                    tibble()
                )
            },
            simplify=TRUE,
            USE.NAMES=FALSE
        ) 
    ) %>% 
    # Now group samples by condition, 
    # these are the groups being comapred to find differential cotacts
    mutate(
        Sample.Group=glue('{Edit}.{Genotype}.{Celltype}'),
        Sample.Group.copy=Sample.Group
    ) %>% 
    group_by(
        isMerged,
        Sample.Group.copy
    ) %>%
    nest(
        samples.df=
            c(
                Sample.Group,
                filepath,
                SampleID,
                resolution,
                Edit,
                Genotype,
                CloneID,
                TechRepID,
                Celltype
            )
    ) %>% 
    ungroup() %>% 
    rename('Sample.Group'=Sample.Group.copy) %>% 
    # Now get all possible pairs of sample groups (project + genotype + CellID
    get_all_row_combinations(
        .,
        {.},
        cols_to_pair=c('isMerged'),
        keep_self=FALSE
    ) %>% 
    # Now each row represents a single compairson of 2 sample groups with all samples 
    mutate(
        samples.df=
            pmap(
                .l=list(samples.df.A, samples.df.B),
                bind_rows
            )
    ) %>%
    select(-c(starts_with('samples.df.'))) %>% 
    # Only compute results with the max resolution of all usable resolutions for individual samples
    rowwise() %>% 
    mutate(
        resolution.max=max(samples.df$resolution),
        resolution.min=min(samples.df$resolution),
    ) %>% 
    ungroup() %>% 
    pivot_longer(
        starts_with('resolution.'),
        names_to='resolution.type',
        names_prefix='resolution.',
        values_to='resolution'
    ) %>% 
    # list all chromosomes separately,
    join_all_rows(tibble(chr=CHROMOSOMES))
}

sample_group_priority_fnc_16p <- function(Sample.Group){
    case_when(
        grepl('16p.DUP.NSC', Sample.Group) ~ 1,  # always numerator in FCs
        grepl('16p.DUP.iN',  Sample.Group) ~ 2,
        grepl('16p.DEL.NSC', Sample.Group) ~ 3,
        grepl('16p.DEL.iN',  Sample.Group) ~ 4,
        grepl('16p.WT.NSC',  Sample.Group) ~ 5,
        grepl('16p.WT.iN',   Sample.Group) ~ 6,
        TRUE ~ -Inf
    )
}

sample_group_priority_fnc_NIPBLWAPL <- function(Sample.Group){
    case_when(
        grepl('WAPL.DEL',  Sample.Group) ~ 1,
        grepl('NIPBL.DEL', Sample.Group) ~ 2,
        grepl('All.WT',    Sample.Group) ~ 3,
        grepl('WAPL.WT',   Sample.Group) ~ 4,
        grepl('NIPBL.WT',  Sample.Group) ~ 5,
        TRUE ~ -Inf
    )
}

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
    covariates.df=NULL,
    zero.p,
    A.min,
    resolution,
    range1,
    range2,
    remove.regions,
    md_plot_file,
    frac.cutoff=0.8,
    effect.col='Sample.Group',
    p.method='fdr',
    ...){
    # row_index=357; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]; md_plot_file=tmp$md_plot_file[[row_index]]; remove.regions=hg38_cyto; p.method='fdr'; effect.col='Sample.Group'; zero.p=tmp$zero.p[[row_index]]; A.min=tmp$A.min[[row_index]]
    # tmp[row_index,]
    # Handle covariates if specified
    design.info <- 
        handle_covariates(
            samples.df,
            covariates.df,
            effect.col
        )
    # print(design.info)
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
        select(-c(resolution)) %>% 
        pmap(
            .l=,
            .f=load_mcool_file,
            resolution=resolution,
            range1=range1,
            range2=range2,
            normalization="NONE",
            .progress=FALSE
        )
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
        # samples.contacts %>% lapply(nrow) %>% unlist()
    # Make experiment object with relevant data+parameters
    make_hicexp(
        data_list=samples.contacts,
        groups=samples.df %>% pull(Sample.Group),
        covariates=design.info$covariates,
        zero.p=zero.p,
        A.min=A.min,
        remove.regions=remove.regions,
        remove_zeroes=FALSE,
        filter=TRUE
    ) %>% 
    # Normalize hic data, use cyclic loess with automatically calculated span
    cyclic_loess(
        verbose=TRUE,
        parallel=TRUE,
        span=NA
    ) %T>% 
    # Plot normalized IFs for all sample pairs
    {
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
    # Preform differential testing on  contacts 
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
    covariates.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    ...){
    # force_redo=TRUE; remove.regions=hg38_cyto; sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL; p.method='fdr'
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    join_all_rows(hyper.params.df) %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    mutate(
        across(
            matches('Sample.Group.(A|B)'),
            ~ sample_group_priority_fnc(.x),
            .names='{.col}.Priority'
        )
    ) %>%
    mutate(
        Sample.Group.Numerator=
            case_when(
                Sample.Group.A.Priority < Sample.Group.B.Priority ~ Sample.Group.A,
                Sample.Group.A.Priority > Sample.Group.B.Priority ~ Sample.Group.B,
                TRUE ~ NA
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == Sample.Group.A ~ Sample.Group.B,
                Sample.Group.Numerator == Sample.Group.B ~ Sample.Group.A,
                TRUE ~ NA
            ),
    ) %>% 
    select(-c(matches('Sample.Group.(A|B)'))) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                MULTIHICCOMPARE_DIR,
                'results',
                glue('merged_{isMerged}'),
                glue('zero.p_{zero.p}'),
                glue('A.min_{A.min}'),
                glue('resolution_{scale_numbers(resolution)}'),
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
        # {.} -> tmp
    # pmap(
    future_pmap(
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
                    ... # arguments taken from columns of .
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load resutls
###################################################
load_multiHiCCompare_results <- function(
    filepath,
    nom.threshold,
    fdr.threshold,
    ...){
    filepath %>%
    read_tsv(show_col_types=FALSE) %>%
    filter(p.adj < fdr.threshold) %>%
    filter(p.value < nom.threshold)
}

load_all_multiHiCCompare_results <- function(
    file_suffix='-multiHiCCompare.tsv',
    sample_group_priority_fnc=NULL,
    resolutions=NULL,
    fdr.threshold=1,
    nom.threshold=1,
    ...){
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
        names=c('Sample.Group.A', 'Sample.Group.B')
    ) %>% 
    mutate(
        across(
            starts_with('Sample.Group'),
            sample_group_priority_fnc,
            .names='{.col}.Priority'
        )
    ) %>% 
    mutate(
        Sample.Group.Numerator=
            case_when(
                Sample.Group.A.Priority < Sample.Group.B.Priority ~ Sample.Group.A,
                Sample.Group.A.Priority > Sample.Group.B.Priority ~ Sample.Group.B
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == Sample.Group.A ~ Sample.Group.B,
                Sample.Group.Numerator == Sample.Group.B ~ Sample.Group.A
            )
    ) %>% 
    select(-ends_with(c('.A', '.B', '.Priority'))) %>% 
    # Only load results with specific params
    mutate(
        results=
            pmap(
                .l=.,
                .f=load_multiHiCCompare_results,
                fdr.threshold=fdr.threshold,
                nom.threshold=nom.threshold,
                .progress=TRUE
            )
    ) %>%
    select(-c(filepath)) %>% 
    unnest(results)
}

###################################################
# Plot Results
###################################################
multiHiCCompare_genome_volcano_plot <- function(
    plot.df,
    color='distance.discrete',
    shape='Comparison',
    scales='fixed',
    pal.direction=-1,
    pal='A',
    size=1,
    alpha=0.6,
    ncol=1,
    ...){
    plot.df %>%
    ggplot() +
    geom_jitter(
        aes(
            x=logFC,
            y=log.p.adj,
            color=.data[[color]],
            shape=.data[[shape]]
            # color=Comparison, size=distance.kb
        ), 
        alpha=alpha
    ) +
    {
        if (is.numeric(plot_df[[color]])) {
            scale_color_viridis(direction=pal.direction, option=pal)
        } else {
            scale_color_viridis(direction=pal.direction, option=pal, discrete=TRUE)
        }
    } +
    facet_wrap(
        ~ Resolution,
        ncol=ncol,
        scales=scales
    ) +
    geom_hline(yintercept=-log10(1e-8), linetype='dashed', color='black') +
    labs(
        title=glue('All Genome Bins with Differential Contacts'),
        color='Contact Distance (Kb)'
    ) +
    scale_y_continuous(
        expand=c(0.01, 0.01, 0.01, 0.01),
        limits=c(0, NA)
    ) +
    theme(
        axis.text.x=element_text(angle=45),
        legend.position='right'
    ) +
    add_ggtheme()
}

multiHiCCompare_allChromosome_volcano_plot <- function(
    plot_df,
    color='distance.discrete',
    shape='Comparison',
    scales='free_y',
    pal='A',
    pal.direction=-1,
    size=1,
    alpha=0.6,
    ncol=1,
    ...){
    g <- 
        plot.df %>%
        ggplot() +
        geom_jitter(
            aes(
                x=logFC,
                y=log.p.adj,
                color=.data[[color]],
                shape=.data[[shape]]
                # color=Comparison, size=distance.kb
            ), 
            alpha=alpha
        ) +
        {
            if (is.numeric(plot_df[[color]])) {
                scale_color_viridis(direction=pal.direction, option=pal)
            } else {
                scale_color_viridis(direction=pal.direction, option=pal, discrete=TRUE)
            }
        } +
        facet_wrap(
            ~ Chr,
            ncol=ncol,
            scales=scales
        ) +
        # geom_hline(yintercept=-log10(1e-8), linetype='dashed', color='black') +
        labs(
            title=glue('All Genome Bins with Differential Contacts'),
            color='Contact Distance (Kb)'
        ) +
        scale_y_continuous(
            expand=c(0.01, 0.01, 0.01, 0.01),
            limits=c(0, NA)
        ) +
        theme(
            strip.text=element_text(size = 20, face='bold'),
            axis.text.x=element_text(angle=45),
            legend.position='right'
        ) +
        add_ggtheme()
    shift_legend(g)
}

multiHiCCompare_chromosome_volcano_plot <- function(
    plot.df,
    chr,
    color='distance.discrete',
    shape='Comparison',
    size=1,
    alpha=0.6,
    scales='fixed',
    ncol=1,
    ...){
    plot.df %>%
    ggplot() +
    geom_jitter(
        aes(
            x=logFC,
            y=log.p.adj,
            color=.data[[color]],
            shape=.data[[shape]]
        ), 
        alpha=alpha,
        size=size
    ) +
    facet_wrap(
        ~ Resolution,
        ncol=ncol,
        scales=scales
    ) +
    labs(
        title=glue('Chr{chr} Bins with Differential Contacts'),
        color='Contact Distance (Kb)'
    ) +
    scale_y_continuous(
        expand=c(0.01, 0.01, 0.01, 0.01),
        limits=c(0, NA)
    ) +
    theme(
        axis.text.x=element_text(angle=45),
        legend.position='right'
    ) +
    add_ggtheme()
}

multiHiCCompare_manhattan_plot <- function(
    plot.df, 
    chr,
    Resolution,
    y_axis='region1.bin',
    n.breaks=50,
    ...){
    plot.df %>%
    ggplot() +
    geom_jitter(
        aes(
            x=value,
            y=.data[[y_axis]],
            color=Comparison
        ), 
    ) +
    facet_wrap(
        ~ statistic,
        nrow=1,
        scales='free_x'
    ) +
    labs(
        title=glue('Chr{chr} Bins with Differential Contacts'),
        x=glue('Chr{chr} bins at {Resolution}  resolution')
    ) +
    scale_y_continuous(
        n.breaks=n.breaks,
        expand=c(0.01, 0.01, 0.01, 0.01),
        limits=c(0, NA)
    ) +
    theme(
        axis.text.x=element_text(angle=45),
        legend.position='right'
    ) +
    add_ggtheme()
}

multiHiCCompare_genome_scatter_plot <- function(
    plot.df,
    size=1,
    alpha=0.6,
    ncol=1,
    ...){
    plot.df %>% 
    ggplot() +
    geom_point(
        aes(
            x=distance.value,
            y=logFC,
            shape=Comparison
        ), 
        alpha=alpha,
        size=size
    ) +
    facet_wrap(
        ~ distance.unit,
        nrow=1,
        scales='free_x'
    ) +
    labs(
        title=glue('Differential Contacts Genome-wide'),
        x='Bin Distance'
    ) +
    scale_y_continuous(
        expand=c(0.01, 0.01, 0.01, 0.01),
        limits=c(0, NA)
    ) +
    theme(
        legend.position='right'
    ) +
    add_ggtheme()
}

