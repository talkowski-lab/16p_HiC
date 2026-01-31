###################################################
# Dependencies
###################################################
# library(tidyverse)
# library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(GGally)
library(scales)

###################################################
# Handling/formatting plots
###################################################
make_ggtheme <- function(...){
    theme(
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            ),
        ...
    )
}

scale_x_axis <- function(
    figure,
    scale.mode='',
    log.base=10,
    axis.label.accuracy=0.1,
    n.breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # Scale y axis based on scale.mode argumetn
    if (scale.mode == 'pct') {
        figure +
        coord_cartesian(xlim=limits) +
        scale_x_continuous(
            expand=expand,
            n.breaks=n.breaks,
            labels=label_percent(),
            ...
        )
    } else if (scale.mode == 'mb') {
        figure +
        coord_cartesian(xlim=limits) +
        scale_x_continuous(
            expand=expand,
            n.breaks=n.breaks,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    } else if (scale.mode == 'log10') {
        figure +
        coord_cartesian(xlim=limits) +
        scale_x_log10(
            expand=expand,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log.base,
                    digits=max(1, -log10(axis.label.accuracy)),
                    signed=FALSE
                ),
            ...
        )
    } else if (scale.mode == 'discrete') {
        figure +
        scale_x_discrete(expand=expand)
    } else if (scale.mode == '') {
        if (is.null(limits)) {
            figure + 
            scale_x_continuous(
                labels=
                    function(x) {
                        format(x, digits=max(1, -log10(axis.label.accuracy)))
                    }
            )
        } else {
            figure + 
            coord_cartesian(xlim=limits) +
            scale_x_continuous(
                expand=expand,
                ...
            )
        }
    }
}

scale_y_axis <- function(
    figure,
    scale.mode='',
    log.base=10,
    axis.label.accuracy=0.1,
    n.breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # Scale y axis based on scale_mode argumetn
    if (scale.mode == 'pct') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            expand=expand,
            n.breaks=n.breaks,
            labels=label_percent(),
            ...
        )
    } else if (scale.mode == 'mb') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            expand=expand,
            n.breaks=n.breaks,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    } else if (scale.mode == 'log10') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_log10(
            expand=expand,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log.base,
                    digits=max(1, -log10(axis.label.accuracy)),
                    signed=FALSE
                ),
            ...
        )
    } else if (scale.mode == 'discrete') {
        figure +
        scale_y_discrete(expand=expand)
    } else if (scale.mode == '') {
        if (is.null(limits)) {
            figure + 
            scale_y_continuous(
                labels=
                    function(x) {
                        format(x, digits=max(1, -log10(axis.label.accuracy)))
                    }
            )
        } else {
            figure + 
            coord_cartesian(ylim=limits) +
            scale_y_continuous(
                expand=expand,
                ...
            )
        }
    }
}

add_faceting <- function(
    figure,
    facet.group=NULL,
    facet.col=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.ncol=NULL,
    ...){
    # Facet as specified
    if (!is.null(facet.col) & !is.null(facet.row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet.row)),
                cols=vars(!!sym(facet.col)),
                ...
            )
    } else if (!is.null(facet.row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet.row)),
                ...
            )
    } else if (!is.null(facet.col)) {
        figure <- 
            figure +
            facet_grid2(
                cols=vars(!!sym(facet.col)),
                ...
            )
    } else if (!is.null(facet.group)) {
        figure <- 
            figure +
            facet_wrap2(
                vars(!!sym(facet.group)),
                nrow=facet.nrow,
                ncol=facet.ncol,
                ...
            )
    }
    figure
}

post_process_plot <- function(
    figure,
    theme.obj=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.col=NULL,
    facet.ncol=NULL,
    facet.group=NULL,
    scales='fixed',
    x.scale.mode='',
    x.log.base=10,
    x.axis.label.accuracy=0.1,
    x.n.breaks=NULL,
    x.limits=NULL,
    x.expand=c(0.00, 0.00, 0.00, 0.00),
    y.scale.mode='',
    y.log.base=10,
    y.axis.label.accuracy=0.1,
    y.n.breaks=NULL,
    y.limits=NULL,
    y.expand=c(0.00, 0.00, 0.00, 0.00),
    plot.elements=NULL,
    ...){
    figure %>% 
    add_faceting(
        scales=scales,
        facet.row=facet.row,
        facet.col=facet.col,
        facet.group=facet.group,
        facet.nrow=facet.nrow,
        facet.ncol=facet.ncol
    ) %>% 
    # Set x-axis scaling (log, Mb, percent etc.)
    scale_x_axis(
        scale.mode=x.scale.mode,
        log.base=x.log.base,
        axis.label.accuracy=x.axis.label.accuracy,
        n.breaks=x.n.breaks,
        limits=x.limits,
        expand=x.expand
    ) %>% 
    # Set y-axis scaling (log, Mb, percent etc.)
    scale_y_axis(
        scale.mode=y.scale.mode,
        log.base=y.log.base,
        axis.label.accuracy=y.axis.label.accuracy,
        n.breaks=y.n.breaks,
        limits=y.limits,
        expand=y.expand
    ) %>% 
    # Add theme elements, as either an object or individual args
    {
        if (!is.null(theme.obj)) {
            . + theme.obj
        } else {
            . + make_ggtheme() + theme(...)
        } 
    } %>% 
    # Add extra elements
    {
        if (!is.null(plot.elements)) {
            reduce(
                .x=plot.elements,
                .f=function(x, y) { x + y },
                .init=.
            )
        } else {
            .
        }
    }
}

###################################################
# Make tabs per plot in Rmd
###################################################
plot_figure_tabs <- function(
    plot.df,
    group.col,
    plot.fnc,
    header.lvl,
    nl.delim,
    return.figure=FALSE,
    merge.base.layers=FALSE,
    grob.nrow=1,
    grob.ncol=NULL,
    ...){
    # List all the groups to plot
    plot.df[[group.col]] %>% 
    as.factor() %>% 
    droplevels() %>% 
    levels() %>%
    # make plot with options for each group of the data
    sapply(
        function(group.value, plot.df, plot.fnc, header.lvl, group.col, nl.delim, return.figure){
            figure <- 
                plot.df %>%
                filter(!!sym(group.col) == group.value) %>%
                plot.fnc(...)
            if (!return.figure & !merge.base.layers) {
                cat(
                    strrep('#', header.lvl), group.value,
                    nl.delim
                )
                print(figure)
                cat(nl.delim)
            } else {
                figure
            }
        },
        plot.df=plot.df,
        plot.fnc=plot.fnc,
        header.lvl=header.lvl,
        group.col=group.col,
        nl.delim=nl.delim,
        return.figure=return.figure,
        simplify=FALSE,
        USE.NAMES=TRUE
    ) %>%
    # Plot in separate tabs unless specified, then make a single, multi-panel plot 
    {
        if (merge.base.layers) {
            cat(
                strrep('#', header.lvl),
                nl.delim
            )
            cowplot::plot.grid(
                plotlist=.,
                nrow=grob.nrow,
                ncol=grob.ncol,
                labels=names(.),
                axis='tb',
                align='hv'
            ) %>%
            print()
            cat(nl.delim)
        } else {
            NULL
        }
    }
}

make_tabs_recursive <- function(
    plot.df, 
    group.cols,
    current.header.lvl,
    plot.fnc,
    tabset.format,
    nl.delim,
    return.figure,
    merge.base.layers,
    ...){
    if (length(group.cols) > 1) {
        group.col <- group.cols[1]
        # message(paste(current.header.lvl, group.col, collapse=','))
        group.values <- 
            plot.df[[group.col]] %>% 
            as.factor() %>% 
            droplevels() %>% 
            levels()
        for (group.value in group.values) {
            cat(
                strrep('#', current.header.lvl), group.value, tabset.format,
                nl.delim
            )
            make_tabs_recursive(
                plot.df=plot.df %>% filter(get({{group.col}}) == group.value),
                group.cols=group.cols[2:length(group.cols)],
                current.header.lvl=current.header.lvl + 1,
                plot.fnc=plot.fnc,
                tabset.format=tabset.format,
                nl.delim=nl.delim,
                return.figure=return.figure,
                merge.base.layers=merge.base.layers,
                ...
            )
        }
    } else if (length(group.cols) == 1) {
        plot_figure_tabs(
            plot.df=plot.df, 
            group.col=group.cols[1],
            header.lvl=current.header.lvl,
            plot.fnc=plot.fnc,
            nl.delim=nl.delim,
            return.figure=return.figure,
            merge.base.layers=merge.base.layers,
            ...
        )
    } else if (length(group.cols) == 0) {
        figure <- 
            plot.df %>%
            plot.fnc(...)
        if (return.figure) {
            return(figure)
        } else {
            cat(
                strrep('#', header.lvl),
                nl.delim
            )
            print(figure)
            cat(nl.delim)
        }
    }
}

make_nested_plot_tabs <- function(
    plot.df,
    group.cols,
    plot.fnc,
    max.header.lvl=2,
    add.top.layer=FALSE,
    tabset.format="{.tabset}",
    nl.delim="\n\n\n",
    return.figure=FALSE,
    merge.base.layers=FALSE,
    ...){
    cat(nl.delim)
    if (add.top.layer) {
        cat(strrep('#', max.header.lvl), tabset.format, nl.delim)
        max.header.lvl <- max.header.lvl + 1
    }
    plot.df %>% 
    make_tabs_recursive(
        group.cols=group.cols,
        current.header.lvl=max.header.lvl,
        plot.fnc=plot.fnc,
        tabset.format=tabset.format,
        nl.delim=nl.delim,
        return.figure=return.figure,
        merge.base.layers=merge.base.layers,
        ...
    )
    cat(nl.delim)
}

###################################################
# Basic Plots
###################################################
plot_barplot <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var='', 
    position='dodge',
    legend.cols=1,
    ...){
    # Set fill group if specified
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    # make it a boxplot 
    { 
        . + 
        geom_col(position=position) +
        guides(fill=guide_legend(ncol=legend.cols))
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_boxplot <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var=NULL, 
    outlier.size=1,
    ...){
    # Set fill group if specified
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    # make it a boxplot 
    { . + geom_boxplot(outlier.size=1) } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_violin <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var=NULL, 
    plot_pts=FALSE,
    jitter.size=0.1,
    quantile.color=NULL,
    draw.quantiles=0L,
    quantile.linewidth=NULL,
    position='dodge',
    adjust=0.5,
    ...){
    # Set fill group if specified
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    # make it a violin plot
    { 
        . + 
        geom_violin(
            quantile.color=quantile.color,
            quantile.linetype=draw.quantiles,
            quantile.linewidth=quantile.linewidth,
            position=position,
            adjust=adjust
        )
    } %>% 
    # plot individual points if specified
    {
        if (plot_pts){
            . +
            geom_jitter(size=jitter.size)
            # geom_jitter(aes(size=jitter.size))
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_heatmap <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var='', 
    label.var=NULL,
    label.size=2,
    label.color='white',
    ...){
    {
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
                fill=.data[[fill.var]]
            )
        ) +
        geom_tile()
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        x.scale.mode='discrete',
        y.scale.mode='discrete',
        ...
    ) %>% 
    {
        if (!(is.null(label.var))) {
            . +
            geom_text(
                aes(label=.data[[label.var]]), 
                color=label.color,
                size=label.size,
            )
        } else {
            .
        }
    }
}

plot_jitter <- function(
    plot.df,
    x.var='',
    y.var='',
    color.var=NULL, 
    shape.var=NULL,
    alpha=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # Set fill group if specified
    {
        if (!is.null(color.var) && !is.null(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    shape=.data[[shape.var]],
                    color=.data[[color.var]]
                )
            )
        } else if (!is.null(color.var) &&  is.null(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    color=.data[[color.var]]
                )
            )
        } else if ( is.null(color.var) && !is.null(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    shape=.data[[shape.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        }
    } %>% 
    # make it a boxplot 
    { 
        . + 
        geom_jitter(
            alpha=alpha,
            size=size
        )
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

###################################################
# gghic utils
###################################################
load_TADs_for_gghic <- function(){
    HITAD_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    post_process_hiTAD_TAD_results() %>% 
    add_column(weight='balanced') %>% 
    dplyr::select(-c(length)) %>% 
    mutate(chr2=chr) %>% 
    nest(TADs=c(chr, start, end)) %>% 
    dplyr::rename(
        'chr'=chr2,
        'TAD.method'=method,
        'SampleID'=Sample.Group
    )
}

load_loops_for_gghic <- function(){
    COOLTOOLS_LOOPS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>%
    post_process_cooltools_dots_results() %>% 
    filter(kernel == 'donut') %>% 
    filter(log10.qval > -log10(0.1)) %>% 
    dplyr::rename(
        'start.P1'=anchor.left,
        'start.P2'=anchor.right
    ) %>% 
    mutate(
        end.P1=start.P1 + resolution,
        end.P2=start.P2 + resolution,
        chr.P1=chr,
        chr.P2=chr,
    ) %>% 
    dplyr::select(
        weight, resolution, SampleID, chr,
        chr.P1, start.P1, end.P1, 
        chr.P2, start.P2, end.P2,
    ) %>% 
    nest(
        loops=
            c(
                chr.P1, start.P1, end.P1, 
                chr.P2, start.P2, end.P2,
            )
    )
}

load_tracks_for_gghic <- function(){
    return(NULL)
}

make_gghic_plot_obj <- function(
    cooler_path, resolution, 
    focus,
    TADs, loops, tracks=NULL,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[1]]', sep='', collapse=';')
    # cooler_path=tmp$cooler_path[[1]];Edit=tmp$Edit[[1]];Celltype=tmp$Celltype[[1]];Genotype=tmp$Genotype[[1]];SampleID=tmp$SampleID[[1]];resolution=tmp$resolution[[1]];focus=tmp$focus[[1]];TAD.method=tmp$TAD.method[[1]];TADs=tmp$TADs[[1]];weight=tmp$weight[[1]];loops=tmp$loops[[1]]
    # make object with contacts
    cc <- 
        ChromatinContacts(
            cooler_path=cooler_path,
            focus=focus,
            resolution=as.integer(resolution)
        ) %>%
        import()
    # add TAD annotations
    if (!is.null(TADs)){
        features(cc, 'TADs') <- GRanges(TADs)
    }
    # add loop annotations
    if (!is.null(loops)){
        features(cc, 'loops') <- 
            GInteractions(  
                anchor1=
                    loops %>% 
                    dplyr::select(ends_with('.P1')) %>% 
                    rename_with(~ str_remove(.x, '.P1')) %>% 
                    GRanges(),
                anchor2=
                    loops %>% 
                    dplyr::select(ends_with('.P2')) %>% 
                    rename_with(~ str_remove(.x, '.P2')) %>% 
                    GRanges()
            ) %>% 
            pairs() %>% 
            makeGInteractionsFromGRangesPairs()
    }
    # add 1D track data
    if (!is.null(tracks)){
        features(cc, 'tracks') <- GRanges(tracks)
    } 
    # return annotated contacts plot object
    cc
}

set_up_gghic_plot_param_sets <- function(
    regions.df,
    hyper.params.df,
    tads.df=NULL,
    loops.df=NULL,
    tracks.df=NULL,
    ...){
    # tads.df=load_TADs_for_gghic(); loops.df=load_loops_for_gghic(); tracks.df=load_tracks_for_gghic(); 
    # list all contact matrices
    list_mcool_files() %>%
    filter(isMerged) %>% 
    dplyr::select(-c(CloneID, TechRepID,  ReadFilter, isMerged)) %>% 
    dplyr::rename(cooler_path=filepath) %>% 
    mutate(SampleID=str_remove(SampleID, '.Merged.Merged')) %>% 
    cross_join(tibble(chr=unique(regions.df$chr))) %>% 
    cross_join(hyper.params.df) %>% 
    # Join TADs
    {
        if (!is.null(tads.df)) {
            left_join(
                .,
                tads.df,
                by=
                    join_by(
                        weight,
                        resolution,
                        SampleID,
                        chr
                    )
            )
        } else {
            add_column(., tads.df=NULL)
        }
    } %>% 
    # Join Loops
    {
        if (!is.null(loops.df)) {
            left_join(
                .,
                loops.df,
                by=
                    join_by(
                        weight,
                        resolution,
                        SampleID,
                        chr,
                    )
            )
        } else {
            add_column(., loops.df=NULL)
        }
    } %>% 
    # Join tracks
    {
        if (!is.null(tracks.df)) {
            left_join(
                .,
                tracks.df,
                by=
                    join_by(
                        resolution,
                        SampleID,
                        chr,
                    )
            )
        } else {
            add_column(., tracks=list(NULL))
        }
    } %>% 
    # Specify regions to plot
    left_join(
        regions.df,
        by=join_by(chr)
    ) %>% 
    # Set up plot params
    rename('resolution.int'=resolution) %>% 
    mutate(resolution=scale_numbers(resolution.int)) %>% 
    mutate(panel.title=glue('{SampleID} @ {resolution}')) %>% 
    arrange(weight, resolution, focus) %>% 
    # Plot all genotypes for each condition together in a single plot
    nest(
        plot.df=
            c(
                Edit, Celltype, Genotype, SampleID,
                cooler_path,
                resolution.int,
                TADs,
                loops,
                tracks,
                panel.title
            )
    ) %>% 
    # need to save as pdf for performance reasons
    mutate(
        # region.title=glue('{SampleID} @ {resolution} in {region.title} ({focus})'),
        plot.title=glue('{region.title} ({focus})'),
        results_file=
            file.path(
                PLOT_DIR,
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('region_{region.title}-gghic.plot.pdf')
            )
    )
}

plot_gghic <- function(
    cooler_path,
    resolution.int,
    focus,
    weight,
    TADs,
    loops,
    tracks,
    panel.title,
    gtf_path=NULL,
    ...){
    # paste(colnames(plots.df), '=plots.df$', colnames(plots.df), '[[1]]', sep='', collapse=';')
    # cooler_path=plots.df$cooler_path[[1]];Edit=plots.df$Edit[[1]];Celltype=plots.df$Celltype[[1]];Genotype=plots.df$Genotype[[1]];SampleID=plots.df$SampleID[[1]];chr=plots.df$chr[[1]];resolution.int=plots.df$resolution.int[[1]];weight=plots.df$weight[[1]];TAD.method=plots.df$TAD.method[[1]];TADs=plots.df$TADs[[1]];loops=plots.df$loops[[1]];tracks=plots.df$tracks[[1]];region.title=plots.df$region.title[[1]];focus=plots.df$focus[[1]];resolution=plots.df$resolution[[1]]
    # make plot obj with contacts + annotation data
    make_gghic_plot_obj(
        cooler_path=cooler_path,
        resolution=resolution.int, 
        focus=focus,
        TADs=TADs, 
        loops=loops,
        tracks=tracks
    ) %>% 
    # Make gghic plot with specified annotations
    {
        gghic(
            .,
            scale_column=weight,
            genome="hg38", 
            gtf_path=gtf_path,
            annotation=!is.null(gtf_path),
            ideogram=TRUE, 
            loop=TRUE,
            tad=TRUE,
            tad_is_0based=TRUE,
            track=TRUE,
            ...
        ) +
        ggtitle(panel.title) +
        theme_hic()
    } #%>% 
    # Handle faceting + scaling + theme options
    # post_process_plot(...)
}

plot_all_regions_gghic <- function(
    plot.df,
    focus,
    weight,
    plot.title,
    output_file=NA,
    width=8,
    height=14,
    ...){
    plot.df %>%
    mutate(
        plot_obj=
            pmap(
                .l=.,
                .f=plot_gghic,
                focus=focus,
                weight=weight,
                ...,
                .progress=TRUE
            )
    ) %>%
    pull(plot_obj) %>% 
    # pmap(
    #     .f=post_process_plot,
    #     ...,
    #     .progress=FALSE
    # ) %>% 
    {
        wrap_plots(
            .,
            ncol=1
        ) +
        plot_annotation(title=unique(plot.df$plot.title))
    } %>%
    {
        if (!is.na(output_file)){
            ggsave(
                output_file, 
                plot=.,
                width=width,
                height=height,
                units='in'
            )
        } else {
            .
        }
    }
}

