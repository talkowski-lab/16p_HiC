###################################################
# Dependencies
###################################################
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
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
    jitter.size=1,
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
            geom_jitter(aes(size=jitter.size))
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
    alpha=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # Set fill group if specified
    {
        if (is.null(color.var)) {
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
                    color=.data[[color.var]]
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

plot_upset <- function(
    plot.df,
    make.binary=FALSE,
    category_col='Comparison',
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

###################################################
# Contact Heatmaps
###################################################
format_plot_params <- function(
    region.df,
    plot_dir,
    title.prefix='RGD Region',
    ...){
    load_mcool_files(
        return_metadata_only=TRUE,
        pattern='*.mapq_30.1000.mcool',
        region.df=region.df,
        range1s=NULL,
        range2s=NULL,
        progress=TRUE
    ) %>% 
    # format querys for loading contacts via fetch()
    rowwise() %>% 
    mutate(
        range1=glue('{region.chr}:{max(0, region.start - window.size)}-{region.end + window.size}'),
        range2=range1
    ) %>%
    ungroup() %>% 
    mutate(
        cis=TRUE,
        region.title=glue('{title.prefix} {region} {region.UCSC}'),
        output_dir=
            file.path(
                PLOT_DIR, 
                glue('region_{region}'),
                glue('normalization_{normalization}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('context_{scale_numbers(window.size)}')
            )
    )
    # group_by(Sample.ID, region) %>% 
    # add_tally(wt=IF, name='region.total.contacts') %>% 
    # ungroup() %>%
    # group_by(Sample.ID, region) %>% 
    # add_count(name='region.nbins') %>% 
    # ungroup() %>%
}

make_contact_plot_df <- function(
    make_diagonal=FALSE,
    make_symmetric=TRUE,
    add_NAs=FALSE,
    fill_var,
    x_var,
    y_var,
    ...){
    load_mcool_file(...) %>% 
    # Explicitly add symetric pairs i.e. 1,2 -> 2,1
    {
        if (make_symmetric) {
            bind_rows(
                .,
                {.} %>% 
                rename_with(
                    .fn=~ case_when(
                        .x == x_var ~ y_var,
                        .x == y_var ~ x_var,
                        TRUE ~ .x
                      )
                ) %>%
                filter(!!sym(x_var) != !!sym(y_var))
            )
        } else {
            . 
        }
    } %>% 
    # Make empty bin-pairs explicit NAs
    {
        if (add_NAs) {
            fill_list <- list(NA)
            names(fill_list) <- fill_var
            group_by(
                .,
                across(
                    -c(
                        x_var,
                        y_var,
                        fill_var
                    )
                )
            ) %>% 
            complete(
                !!sym(x_var),
                !!sym(y_var),
                fill=fill_list
            ) %>% 
            ungroup()
        } else {
            . 
        }
    } %>% 
    # Transoforming for triangle version, stolen from HiContacts library
    # https://github.com/js2264/HiContacts/blob/146485f3e517429f7d4aa408898f6e89d4a04b36/R/plotting.R#L422
    {
        if (make_diagonal) {
            mutate(
                .,
                distance=range2 - range1,
                bin=range1 + (range2 - range1) / 2
            )
        } else {
            .
        }
    }
}

plot_contacts_heatmap <- function(
    contacts,
    resolution,
    x_var='range1',
    y_var='range2',
    fill_var='IF',
    transform_fnc=log10,
    region.start=NULL, 
    region.end=NULL,
    facet_col=NULL,
    facet_row=NULL,
    scales='fixed',
    cmap=coolerColors(),
    axis.label.accuracy=0.01,
    x_text_angle=25,
    linetype='solid',
    linecolor='black',
    xlinewidth=0.7,
    ylinewidth=0.7,
    na.color="white",
    ...){
    figure <- 
        contacts %>%
        ggplot(
            aes(
                x=.data[[x_var]],
                y=.data[[y_var]]
            )
        ) +
        geom_tile(
            aes(fill=transform_fnc(.data[[fill_var]])),
            width=resolution
        ) +
        # geom_tile(aes(fill=fill_data)) +
        scale_fill_gradientn(
            colors=cmap,
            na.value=na.color
        ) +
        scale_x_continuous(
            expand=c(0,0,0,0),
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                )
        ) +
        scale_y_continuous(
            expand=c(0,0,0,0),
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                )
        ) +
        theme(
            axis.text.x=element_text(angle=x_text_angle, hjust=1),
            legend.position='top'
        ) 
    # add vertical lines
    if (xlinewidth > 0) {
        figure <- 
            figure +
            geom_vline(
                xintercept=c(region.start, region.end),
                linetype=linetype,
                color=linecolor,
                linewidth=xlinewidth
            )
    }
    # add horizontal lines
    if (ylinewidth > 0) {
        figure <- 
            figure +
            geom_hline(
                yintercept=c(region.start, region.end),
                linetype=linetype,
                color=linecolor,
                linewidth=ylinewidth
            )
    }
    # add faceting
    figure <- 
        figure %>% 
        add_faceting(
            facet_col=facet_col,
            facet_row=facet_row,
            scales='fixed'
        )
    figure
}

heatmap_wrapper <- function(
    filepath, normalization, resolution, cis, range1, range2,
    make_diagonal=FALSE, make_symmetric=TRUE, add_NAs=FALSE, 
    x_var, y_var, fill_var,
    Sample.ID,
    region.title,
    window.size,
    xlab, ylab, fill_lab,
    output_file, width=7, height=7,
    ...){
    contacts <- 
        make_contact_plot_df(
            make_diagonal=make_diagonal,
            make_symmetric=make_symmetric,
            add_NAs=add_NAs,
            x_var=x_var,
            y_var=y_var,
            fill_var=fill_var,
            filepath=filepath,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2,
            cis=cis
        )
    # Now we can plot the figure
    if (make_diagonal) {
        x_var <- 'bin'
        y_var <- 'distance'
    }
    figure <- 
        contacts %>% 
        plot_contacts_heatmap(
            resolution=resolution,
            x_var=x_var,
            y_var=y_var,
            ...
        ) +
        labs(
            title=glue('{region.title} +/- {scale_numbers(window.size)}'),
            subtitle=glue('{Sample.ID} | {scale_numbers(resolution)} | normalization={normalization}'),
            fill=fill_lab,
            x=xlab,
            y=ylab
        ) +
        make_ggtheme() +
        theme(
            legend.position='right',
            strip.text=element_text(face='bold', size=10),
            axis.title=element_text(face='bold', size=10)
        ) 
    ggsave(
        filename=output_file,
        plot=figure,
        width=width, 
        height=height,
        units='in',
        create.dir=TRUE
    )
}

###################################################
# log2FC Heatmaps
###################################################
make_sample_pairs <- function(
    annotated.contacts.df,
    ...){
    annotated.contacts.df %>%
    nest(data=-c(filepath, Sample.ID)) %>% 
    # Get all pairs of samples with comparable contact matrices
    full_join(
        .,
        {.},
        suffix=c('.A', '.B'),
        by=join_by(data)  # all other params must be equal
    ) %>%
    # Dedup pairs (symetrical)
    filter(Sample.ID.A != Sample.ID.B) %>% 
    rowwise() %>% 
    mutate(pair.id=paste0(sort(c(Sample.ID.A, Sample.ID.B)), collapse="")) %>% 
    ungroup() %>% 
    distinct(
        pair.id, 
        .keep_all=TRUE
    ) %>%
    select(-c(pair.id)) %>% 
    unnest(data)
}

NIPBL_WAPL_sample_priority_fnc <- function(Sample.ID){
    case_when(
        grepl( 'WAPL.DEL', Sample.ID) ~ 1,
        grepl('NIPBL.DEL', Sample.ID) ~ 2,
        grepl( 'WAPL.WT',  Sample.ID) ~ 3,
        grepl('NIPBL.WT',  Sample.ID) ~ 4,
        TRUE ~ -Inf
    )
}

logfc_heatmap_wrapper <- function(
    Sample.ID.A, Sample.ID.B,
    filepath.A, filepath.B,
    sample_priority_fnc,
    contact.params, window.size, resolution, normalization,
    region.title,
    xlab, ylab, fill_lab,
    output_file, width=7, height=7,
    ...){
    # Determine which sample should be numerator
    if (sample_priority_fnc(Sample.ID.A) < sample_priority_fnc(Sample.ID.B)) {
        tmp.ID.holder <- Sample.ID.A
        tmp.filepath <- filepath.A
        Sample.ID.A <- Sample.ID.B
        filepath.A <- filepath.B
        Sample.ID.B <- tmp.ID.holder
        filepath.B <- tmp.filepath
    }
    # Contacts for sample A
    contacts.A <- 
        contact.params %>% 
        pmap(
            .f=make_contact_plot_df,
            filepath=filepath.A
        ) %>%
        first()
    # Contacts for sample B
    contacts.B <- 
        contact.params %>% 
        pmap(
            .f=make_contact_plot_df,
            filepath=filepath.B
        ) %>%
        first()
    # Join and calculate logFC of contacts with appropriate numerator
    contacts <- 
        inner_join(
            contacts.A,
            contacts.B,
            suffix=c('.A', '.B'),
            by=
                join_by(
                    chr,
                    range1,
                    range2
                )
        ) %>%
        mutate(log2fc=log2(IF.A / IF.B))
    # Now we can plot the figure
    figure <- 
        contacts %>% 
        plot_contacts_heatmap(
            resolution=resolution,
            fill_var='log2fc',
            ...
        ) +
        labs(
            title=glue('{Sample.ID.A} vs {Sample.ID.B}'),
            caption=glue('{region.title} +/- {scale_numbers(window.size)}\nresolution={scale_numbers(resolution)}\nnormalization={normalization}'),
            fill=fill_lab,
            x=xlab,
            y=ylab
        ) +
        make_ggtheme() +
        theme(
            legend.position='right',
            strip.text=element_text(face='bold', size=10),
            axis.title=element_text(face='bold', size=10)
        ) 
    ggsave(
        filename=output_file,
        plot=figure,
        width=width, 
        height=height,
        units='in',
        create.dir=TRUE
    )
}

