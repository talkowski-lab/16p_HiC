# Dependencies
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(scales)
###############
# Handling/formatting plots
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

scale_y_axis <- function(
    figure,
    mode='',
    log_base=10,
    axis_label_accuracy=2,
    n_breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # Scale y axis based on mode argumetn
    if (mode == 'pct') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            # limits=limits,
            expand=expand,
            n.breaks=n_breaks,
            labels=label_percent(),
            ...
        )
    } else if (mode == 'mb') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            # limits=limits,
            expand=expand,
            n.breaks=n_breaks,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis_label_accuracy
                ),
            ...
        )
    } else if (mode == 'log10') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_log10(
            # limits=limits,
            expand=expand,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log_base,
                    digits=axis_label_accuracy,
                    signed=FALSE
                ),
            ...
        )
    } else if (mode == '') {
        if (is.null(limits)) {
            figure + 
            scale_y_continuous(
                labels=
                    function(x) {
                        format(x, digits=axis_label_accuracy)
                    }
            )
        } else {
            figure + 
            coord_cartesian(ylim=limits) +
            scale_y_continuous(
                # limits=limits,
                expand=expand,
                ...
            )
        }
        # coord_cartesian(ylim=limits) +
        # scale_y_continuous(...)
    }
}

add_faceting <- function(
    figure,
    facet_group=NULL,
    facet_col=NULL,
    facet_row=NULL,
    ...){
    # Facet as specified
    if (!is.null(facet_col) & !is.null(facet_row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                cols=vars(!!sym(facet_col)),
                ...
            )
    } else if (!is.null(facet_row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                ...
            )
    } else if (!is.null(facet_col)) {
        figure <- 
            figure +
            facet_grid2(
                cols=vars(!!sym(facet_col)),
                ...
            )
    } else if (!is.null(facet_group)) {
        figure <- 
            figure +
            facet_wrap2(
                vars(!!sym(facet_group)),
                ...
            )
    }
    figure
}

post_process_plot <- function(
    figure,
    theme_obj=NULL,
    # add_theme=TRUE,
    facet_row=NULL,
    facet_col=NULL,
    facet_group=NULL,
    scales='fixed',
    scale_mode='',
    log_base=10,
    axis_label_accuracy=2,
    n_breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    # legend.position='right',
    # legend.ncols=1,
    # axis.text.x=element_text(angle=35, hjust=1),
    ...){
    figure <- 
        add_faceting(
            figure,
            facet_row=facet_row,
            facet_col=facet_col,
            facet_group=facet_group,
            scales=scales
        )
    # Set y-axis scaling (log, Mb, percent etc.)
    figure <- 
        figure %>% 
        scale_y_axis(
            mode=scale_mode,
            log_base=log_base,
            axis_label_accuracy=axis_label_accuracy,
            n_breaks=n_breaks,
            limits=limits,
            expand=expand
        )
    # Add theme elements, as either an object or individual args
    if (!is.null(theme_obj)) {
        figure <- figure + theme_obj
    } else {
        figure <- figure + make_ggtheme(...)
    } 
    figure
}
###############
# Make tabs per plot in Rmd
plot_figure_tabs <- function(
    plot.df,
    group_col,
    plot_fnc,
    header_lvl,
    nl_delim,
    return_figure=FALSE,
    ...){
    # message(paste(header_lvl, group_col, collapse=','))
    # print(table(plot.df[[group_col]]))
    plot.df[[group_col]] %>% 
    as.factor() %>% 
    droplevels() %>% 
    levels() %>%
    sapply(
        function(group_value, plot.df, plot_fnc, header_lvl, group_col, nl_delim, return_figure){
            figure <- 
                plot.df %>%
                # filter(get({{group_col}}) == group_value) %>%
                filter(!!sym(group_col) == group_value) %>%
                plot_fnc(...)
                if (return_figure) {
                    figure
                } else {
                cat(
                    strrep('#', header_lvl), group_value,
                    # nl_delim, "Rows per df", nrow(plot.df),
                    nl_delim
                )
                print(figure)
                cat(nl_delim)
            }
        },
        plot.df=plot.df,
        plot_fnc=plot_fnc,
        header_lvl=header_lvl,
        group_col=group_col,
        nl_delim=nl_delim,
        return_figure=return_figure,
        simplify=FALSE,
        USE.NAMES=TRUE
    )
}

make_tabs_recursive <- function(
    plot.df, 
    group_cols,
    current_header_lvl,
    plot_fnc,
    tabset_format,
    nl_delim,
    return_figure=FALSE,
    ...){
    # cat("LENGTH OF GROUP COLS", length(group_cols), group_cols, "\n\n\n")
    if (length(group_cols) == 1) {
        plot_figure_tabs(
            plot.df=plot.df, 
            group_col=group_cols[1],
            header_lvl=current_header_lvl,
            plot_fnc=plot_fnc,
            nl_delim=nl_delim,
            return_figure=return_figure,
            ...
        )
    } else {
        group_col <- group_cols[1]
        # message(paste(current_header_lvl, group_col, collapse=','))
        group_values <- 
            plot.df[[group_col]] %>% 
            as.factor() %>% 
            droplevels() %>% 
            levels()
        for (group_value in group_values) {
            cat(
                strrep('#', current_header_lvl), group_value, tabset_format,
                # nl_delim, "Rows:", nrow(plot.df), 
                nl_delim
            )
            make_tabs_recursive(
                plot.df=plot.df %>% filter(get({{group_col}}) == group_value),
                group_cols=group_cols[2:length(group_cols)],
                current_header_lvl=current_header_lvl + 1,
                plot_fnc=plot_fnc,
                tabset_format=tabset_format,
                nl_delim=nl_delim,
                ...
            )
        }
    }
}

make_nested_plot_tabs <- function(
    plot.df,
    group_cols,
    plot_fnc,
    max_header_lvl=2,
    tabset_format="{.tabset .tabset-pills}",
    nl_delim="\n\n\n",
    return_figure=FALSE,
    ...){
    cat(nl_delim)
    cat(strrep('#', max_header_lvl), tabset_format, nl_delim)
    plot.df %>% 
    make_tabs_recursive(
        group_cols=group_cols,
        current_header_lvl=max_header_lvl+1,
        plot_fnc=plot_fnc,
        tabset_format=tabset_format,
        nl_delim=nl_delim,
        return_figure=return_figure,
        ...
    )
    cat(nl_delim)
}
###############
# Basic Plots
plot_barplot <- function(
    plot.df,
    x_var='',
    y_var='',
    fill_var='', 
    position='dodge',
    legend_cols=1,
    # facet_group=NULL,
    # facet_col=NULL,
    # facet_row=NULL,
    # scales='fixed',
    # scale_mode='',
    # legend.position='right',
    # axis.text.x=element_text(angle=35, hjust=1),
    ...){
    # Set fill group if specified
    figure <- 
        if (is.null(fill_var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x_var]],
                    y=.data[[y_var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x_var]],
                    y=.data[[y_var]],
                    fill=.data[[fill_var]]
                )
            )
        }
    # make it a boxplot 
    figure <- 
        figure + 
        geom_col(position=position) +
        guides(fill=guide_legend(ncol=legend_cols))
    # Handle faceting + scaling + theme options
    figure <- 
        figure %>% 
        post_process_plot(
            # facet_row=facet_row,
            # facet_col=facet_col,
            # facet_group=facet_group,
            # scales=scales,
            # scale_mode=scale_mode,
            # axis_label_accuracy=axis_label_accuracy,
            # legend.position=legend.position,
            # axis.text.x=axis.text.x,
            ...
        )
    figure
}

plot_boxplot <- function(
    plot.df,
    x_var='',
    y_var='',
    fill_var=NULL, 
    facet_row=NULL,
    facet_col=NULL,
    facet_group=NULL,
    scales='fixed',
    scale_mode='',
    y_accuracy=1,
    # legend.position='top',
    # axis.text.x=element_text(angle=55, hjust=1),
    ...){
    # Set fill group if specified
    figure <- 
        if (is.null(fill_var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x_var]],
                    y=.data[[y_var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x_var]],
                    y=.data[[y_var]],
                    fill=.data[[fill_var]]
                )
            )
        }
    # make it a boxplot 
    figure <- 
        figure + 
        geom_boxplot()
    # Handle faceting + scaling + theme options
    figure <- 
        figure %>% 
        post_process_plot(
            facet_row=facet_row,
            facet_col=facet_col,
            facet_group=facet_group,
            scales=scales,
            mode=scale_mode,
            axis_label_accuracy=y_accuracy,
            # legend.position=legend.position,
            # axis.text.x=axis.text.x,
            ...
        )
    figure
}

plot_heatmap <- function(
    plot.df,
    x_var='',
    y_var='',
    fill_var='', 
    facet_row=NULL,
    facet_col=NULL,
    facet_group=NULL,
    scale_mode='',
    scales='fixed',
    # legend.position='right',
    # axis.text.x=element_text(angle=55, hjust=1),
    ...){
    figure <- 
        ggplot(
            plot.df,
            aes(
                x=.data[[x_var]],
                y=.data[[y_var]],
                fill=.data[[fill_var]]
            )
        ) +
        geom_tile()
    # Handle faceting + scaling + theme options
    figure <- 
        figure %>% 
        post_process_plot(
            facet_row=facet_row,
            facet_col=facet_col,
            facet_group=facet_group,
            scales=scales,
            mode=scale_mode,
            axis_label_accuracy=y_accuracy,
            # legend.position=legend.position,
            # axis.text.x=axis.text.x,
            ...
        )
    figure
}
###############
# Contact Heatmaps
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
    axis_label_accuracy=0.01,
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
                    accuracy=axis_label_accuracy
                )
        ) +
        scale_y_continuous(
            expand=c(0,0,0,0),
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis_label_accuracy
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
###############
# log2FC Heatmaps
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
###############
# DEPRECATED
DEP_plot_coverage_lineplot <- function(
    plot.df,
    y_val,
    y_lab,
    region,
    facet_row=NULL, facet_col=NULL,
    scales='fixed', indpt=NULL,
    hl_start=NULL, hl_end=NULL, hl_color='grey',
    n_ticks=12,
    ...){
    # set tickmarks along x axis 
    breaks <- 
        seq(
            min(plot.df$start),
            max(plot.df$start),
            (max(plot.df$start) - min(plot.df$start)) / n_ticks
        )
    labels <- 
        breaks %>%
        {. / 1000} %>% 
        format(
            scientific=FALSE,
            trim=TRUE,
            digist=1
        ) %>%
        paste0('Kb')
    # Make plot
    g <- 
        plot.df %>%
        ggplot(
            aes(
                x=start,
                y=.data[[y_val]],
                color=Sample.ID
            )
        ) +
        # Draw lines for each sample across the region
        geom_path() +
        scale_x_continuous(
            breaks=breaks,
            labels=labels
        )
    # Shade sub-region if specified
    if (!(is.null(hl_start))) {
        g <- 
            g +
            geom_vline(
                xintercept=
                    c(
                        hl_start,
                        hl_end
                    ), 
                    linetype='dashed',
                    linewidth=0.1
            ) +
            # Shade any specified sub-region of interest
            geom_rect(
                ymin=-Inf,
                ymax=Inf,
                xmin=hl_start,
                xmax=hl_end, 
                fill=hl_color,
                alpha=0.1,
                show.legend=FALSE
            )
    }
    # Handle faceting
    g <- 
        g + 
        facet_grid2(
            cols=ifelse(is.null(facet_row), NULL, facet_row),
            rows=ifelse(is.null(facet_row), NULL, facet_row),
            scales=scales,
            independent=indpt
        ) +
        # Visual options
        labs(
            x='Genomic Position',
            y=y_lab
        ) +
        theme(
            legend.position='top',
            axis.text.x=
                element_text(
                    angle=35,
                    hjust=1
                )
        ) +
        make_ggtheme()
    # return plot object
    g
}
