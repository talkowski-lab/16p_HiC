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
# Load resutls
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
    multiHiCCompare_results_dir,
    comparisons=NULL,
    resolutions=NULL,
    fdr.threshold=1,
    nom.threshold=1,
    ...){
    multiHiCCompare_results_dir %>% 
    list.files(
        full.names=FALSE,
        recursive=TRUE,
        pattern='*.tsv'
    ) %>% 
    tibble(filename=.) %>%
    mutate(filepath=file.path(multiHiCCompare_results_dir, filename)) %>%
    # Get info of how results were produced
    mutate(filename=str_remove(basename(filename), '.tsv$')) %>% 
    separate_wider_delim(
        filename,
        cols_remove=FALSE,
        delim='_',
        names=
            c(
                'Comparison',
                'Chr',
                'Resolution',
                'zero.p',
                'A.min',
                'Is.Filtered'
            )
    ) %>%
    mutate(Resolution=as.integer(Resolution)) %>% 
    # Only load results with specific params
    { 
        if (is.null(resolution)) {
            .
        } else {
            filter(., Resolution %in% c(resolutions))
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
    unnest(results) %>% 
    mutate(
        sig.lvl=
            case_when(
                p.adj < 0.0001 ~ 'FDR.1e-3',
                p.adj < 0.001 ~ 'FDR.1e-3',
                p.adj < 0.01 ~ 'FDR.1e-2',
                p.adj < 0.1 ~ 'FDR.1e-1',
                p.value < 0.00005 ~ 'Nom.5e-5',
                p.value < 0.001 ~ 'Nom.5e-3',
                p.value < 0.05 ~ 'Nom.5e-2',
                TRUE ~ 'NS'
            )
    )
}
###############
# Plot Results
multiHiCCompare_genome_volcano_plot <- function(
    plot_df,
    output_file,
    color='distance.discrete',
    shape='Comparison',
    scales='fixed',
    pal.direction=-1,
    pal='A',
    size=1,
    alpha=0.6,
    ncol=1,
    width=7,
    height=7,
    ...){
    g <- 
        ggplot(plot_df) +
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
    ggsave( 
        filename=output_file,
        plot=g,
        width=width,
        height=height,
        units='in'
    )
}
multiHiCCompare_allChromosome_volcano_plot <- function(
    plot_df,
    output_file,
    color='distance.discrete',
    shape='Comparison',
    scales='free_y',
    pal='A',
    pal.direction=-1,
    size=1,
    alpha=0.6,
    ncol=1,
    width=7,
    height=7,
    ...){
    g <- 
        ggplot(plot_df) +
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
    ggsave( 
        filename=output_file,
        plot=shift_legend(g),
        width=width,
        height=height,
        units='in'
    )
}
multiHiCCompare_chromosome_volcano_plot <- function(
    plot_df,
    chr,
    output_file,
    color='distance.discrete',
    shape='Comparison',
    size=1,
    alpha=0.6,
    scales='fixed',
    ncol=1,
    width=7,
    height=7,
    ...){
    g <- 
        ggplot(plot_df) +
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
    ggsave( 
        filename=output_file,
        plot=g,
        width=width,
        height=height,
        units='in'
    )
}
multiHiCCompare_manhattan_plot <- function(
    plot_df, 
    chr,
    Resolution,
    y_axis='region1.bin',
    n.breaks=50,
    width=7,
    height=15,
    output_file,
    ...){
    g <- 
        ggplot(plot_df) +
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
    ggsave( 
        filename=output_file,
        plot=g,
        width=width,
        height=height,
        units='in'
    )
}
multiHiCCompare_genome_scatter_plot <- function(
    plot_df,
    size=1,
    alpha=0.6,
    ncol=1,
    width=7,
    height=7,
    output_file,
    ...){
    g <- 
        ggplot(plot_df) +
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
    ggsave( 
        filename=output_file,
        plot=g,
        width=width,
        height=height,
        units='in'
    )
}
