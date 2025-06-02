# Dependencies
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(scales)
###############
# Formatting
add_ggtheme <- function(){
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
            )
    )
}

shift_legend <- function(p){
    if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
    } else {
    gp <- p
    }

# check for unfilled facet panels
    facet.panels <- grep("^panel", gp[["layout"]][["name"]])
    empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
    empty.facet.panels <- facet.panels[empty.facet.panels]
    if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
    }

# establish extent of unfilled facet panels (including any axis cells in between)
    empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
    empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
    names(empty.facet.panels) <- c("t", "l", "b", "r")

# extract legend & copy over to location of unfilled facet panels
    guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
    if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
    }
    gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")

# squash the original guide box's row / column (whichever applicable)
# & empty its cell
    guide.grob <- gp[["layout"]][guide.grob, ]
    if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
    }
    if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
    }
    gp <- gtable_remove_grobs(gp, "guide-box")

    return(gp)
}

###############
# Utility
make_tab_per_group <- function(
    plot.df, 
    group_col,
    plot_fnc,
    col_header=TRUE,
    lvl=3,
    ...){

    cat("\n\n")
    max_header <- rep('#', max_lvl) %>% paste0(collapse="")
    cat(max_header, ifelse(col_header, group_col, ""), "{.tabset}", "\n\n")
    tab_header <- paste0(max_header, '#')
    group_values <- plot.df[[group_col]] %>% as.factor() %>% droplevels() %>% levels()
    for (group_value in group_values){
        cat(tab_header, group_value, "\n\n")
        print("PLOT GOES HERE")
        figure <- 
            plot.df %>%
            filter(get({{group_col}}) == group_value) %>%
            plot_fnc(...)
        print(figure)
        cat("\n\n")
    }
}

plot_figure_tabs <- function(
    plot.df,
    group_col,
    plot_fnc,
    header_lvl,
    nl_delim,
    ...){
    plot.df[[group_col]] %>% 
        as.factor() %>% 
        droplevels() %>% 
        levels() %>%
        lapply(
            function(group_value, plot.df, plot_fnc, header_lvl, group_col, nl_delim){
                cat(
                    strrep('#', header_lvl), group_value,
                    # nl_delim, "Rows per df", nrow(plot.df),
                    nl_delim
                )
                figure <- 
                    plot.df %>%
                    filter(get({{group_col}}) == group_value) %>%
                    plot_fnc(...)
                print(figure)
                cat(nl_delim)
            },
            plot.df=plot.df,
            plot_fnc=plot_fnc,
            header_lvl=header_lvl,
            group_col=group_col,
            nl_delim=nl_delim
        )
}

make_tabs_recursive <- function(
    plot.df, 
    group_cols,
    current_header_lvl,
    plot_fnc,
    tabset_format,
    nl_delim,
    ...){
    # cat("LENGTH OF GROUP COLS", length(group_cols), group_cols, "\n\n\n")
    if (length(group_cols) == 1) {
        plot_figure_tabs(
            plot.df=plot.df, 
            group_col=group_cols[1],
            header_lvl=current_header_lvl,
            plot_fnc=plot_fnc,
            nl_delim=nl_delim,
            ...
        )
    } else {
        group_col <- group_cols[1]
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
        ...
    )
    cat(nl_delim)
}
###############
# Common Plots
plot_coverage_lineplot <- function(
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
        add_ggtheme()
    # return plot object
    g
}

plot_contacts_heatmap <- function(
    contacts,
    x_var='range1',
    y_var='range2',
    fill_var='IF',
    transform_fnc=log10,
    facet_col=NULL,
    facet_row=NULL,
    scales='fixed',
    cmap=coolerColors(),
    axis_label_accuracy=0.01,
    x_text_angle=25,
    xlinewidth=0.3,
    ylinewidth=0.3,
    cfr=NULL,
    na.value="white",
    ...){
    figure <- 
        contacts %>%
        ggplot(
            aes(
                x=.data[[x_var]],
                y=.data[[y_var]]
            )
        ) +
        geom_tile(aes(fill=transform_fnc(.data[[fill_var]]))) +
        # geom_tile(aes(fill=fill_data)) +
        scale_fill_gradientn(
            colors=cmap,
            na.value=na.value
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
        theme(axis.text.x=element_text(angle=x_text_angle, hjust=1)) 
    # add vertical lines
    if (xlinewidth > 0) {
        figure <- 
            figure +
            geom_vline(
                aes(xintercept=region.start),
                linetype='solid',
                color='black',
                linewidth=xlinewidth
            ) +
            geom_vline(
                aes(xintercept=region.end),
                linetype='solid',
                color='black',
                linewidth=xlinewidth
            )
    }
    # add horizontal lines
    if (ylinewidth > 0) {
        figure <- 
            figure +
            geom_hline(
                aes(yintercept=region.start),
                linetype='solid',
                color='black',
                linewidth=ylinewidth
            ) +
            geom_hline(
                aes(yintercept=region.end),
                linetype='solid',
                color='black',
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
    # change scaling
    if (!(is.null(cfr))) {
        figure <- figure + coord_fixed(ratio=cfr)
    }
    figure
}

heatmap_wrapper <- function(
    contacts,
    Sample.ID,
    region.title,
    resolution,
    window.size,
    xlab,
    ylab,
    fill_lab,
    output_file,
    width=7,
    height=7,
    ...){
    figure <- 
        contacts %>% 
        plot_contacts_heatmap(...) +
        labs(
            title=glue('{region.title} +/- {window.size}'),
            subtitle=glue('{Sample.ID} @{resolution}'),
            fill=fill_lab,
            x=xlab,
            y=ylab
        ) +
        add_ggtheme() +
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
