# Depdendencies
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
###############
# Generate resutls
run_multiHiCCompare <- function(
    sample_df,
    sample_groups,
    region,
    resolution,
    zero.p,
    A.min,
    filter_bins,
    remove.regions=hg38_cyto,
    output_file=NA,
    md_plot_file,
    ...){
    # sample_groups=c('16p11.WT', '16p11.DEL'); region='chr2'; resolution=100000; zero.p=0.8; A.min=5; filter_bins=TRUE; remove.regions=hg38_cyto; output_file=glue('{MULTIHICCOMPARE_DIR}/results/test_results.tsv'); md_plot_file=glue('{MULTIHICCOMPARE_DIR}/MD_plots/test_mdplot.pdf')
    if (!is.na(output_file)) {
        if (file.exists(output_file)) {
            print(glue('cached results exists: {output_file}, skipping'))
            return(NA)
        }
    }
    sample_df %>% 
    dplyr::filter(sample_group %in% sample_groups) %>% 
    # extract contact matrices, each row is a pair of bins + raw interaction freq.
    rowwise() %>% 
    mutate(
        filepath=
            file.path(
                MATRIX_DIR, 
                glue('{sample_ID}.{format({resolution}, scientific=FALSE)}.{region}.txt')
            )
    ) %>%
    mutate(
        sparse.matrix=
            filepath %>% 
            read.table(header=FALSE) %>%
            {.[,c(1,2,5,7)]} %>%
            setNames(c('chr', 'region1', 'region2', 'IF')) %>% 
            mutate(chr=str_remove(chr, '[Cc]hr')) %>% 
            list()
    ) %>%
    ungroup() %>% 
    make_hicexp(
        data_list=.$sparse.matrix,
        groups=.$sample_group,
        filter=filter_bins,
        zero.p=zero.p, # remove rows with (% interactions == 0) >=  zero.p
        A.min=A.min,   # remove rows where (avg interaction freq) < A.min
        remove.regions=remove.regions
    ) %>% 
    # Normalize hic data, use cyclic loess with automatically calculated span
    fastlo(
        verbose=TRUE,
        parallel=TRUE,
        span=NA
    ) %T>% 
    # Plot normalized IFs 
    {
        pdf(md_plot_file, height=nrow(.@metadata) * 2, width=6)
        MD_hicexp(., pcol=2)
        dev.off()
    } %>%
    # Now preform differential testing
    hic_exactTest(
        p.method='fdr',
        parallel=TRUE
    ) %>% 
    results() %>%
    as_tibble() %>%
    arrange(p.adj) %>%
    {
        if (is.na(output_file)){
            .
        } else {
            write_tsv(., output_file)
        }
    }
}

###############
# Load resutls
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
    comparisons=NULL,
    resolutions=NULL,
    fdr.threshold=1,
    nom.threshold=1,
    ...){
    parse_results_filelist(
        input_dir=MULTIHICCOMPARE_DIR,
        suffix='-multiHiCCompare.tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    # Only load results with specific params
    { 
        if (is.null(resolution)) {
            .
        } else {
            filter(., resolution %in% c(resolutions))
        }
    } %>%
    { 
        if (is.null(comparisons)) {
            .
        } else {
            filter(., grepl(comparisons, Comparison))
        }
    } %>%
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
###############
# Plot Results
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
