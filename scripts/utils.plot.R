library(ggplot2)
library(ggpubr)
library(ggh4x)
###############
# Misc
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
# Plot Matrix QC stuff
plot_bin_totals_stats  <- function(plot.df){
    plot.df %>%
    ggplot(
        aes(
            x=SampleID,
            y=total,
            color=SampleID
        )
    ) +
    geom_jitter(size=2) +
    geom_boxplot(oulier=NA) +
    facet_wrap2(~ Chr, scales='fixed') +
    labs(title=glue('{plot.df$normalization[[1]]} normalized count @ {plot.df$resolution[[1]] / 1000}Kb')) +
    add_ggtheme()
}
plot_distance_freqs <- function(
    plot.df,
    y_val='value',
    ...){
    plot.df %>%
    ggplot(aes(x=Category, y=.data[[y_val]], fill=SampleID)) +
    geom_col(position='dodge') +
    scale_y_continuous(breaks=seq(0, max(plot.df[[y_val]]), 5)) +
    # scale_y_continuous(breaks=~ seq(.x, .y, 5)) +
    labs(y='% Unique Read Pairs') +
    theme(axis.text.x=element_text(hjust=1, angle=35)) +
    add_ggtheme()
}
plot_totals_across_chr16 <- function(
    plot.df,
    ...){
    elise.start=29488679; elise.end=30188679
    breaks=
        seq(
            min(plot.df$bin.start),
            max(plot.df$bin.start),
            (max(plot.df$bin.start) - min(plot.df$bin.start)) / 12
        )
    labels=glue('{format(breaks / 1000, scientific=FALSE, trim=TRUE, digits=1)}Kb')
    plot.df %>%
    ggplot(aes(x=bin.start, color=SampleID)) +
    geom_path(aes(y=total)) +
    facet_grid2(cols=vars(normalization), rows=vars(resolution), scales='free_y', independent='y') +
    # geom_rect(ymin=-Inf, ymax=Inf, xmin=elise.start, xmax=elise.end, fill='grey', alpha=0.1, show.legend=FALSE) +
    geom_vline(xintercept=c(elise.start, elise.end), linetype='dashed', linewidth=0.1) +
    scale_x_continuous(breaks=breaks, labels=labels) +
    labs(x='Chromosome bin') +
    theme(legend.position='top', axis.text.x=element_text(angle=35, hjust=1)) +
    add_ggtheme()
}
###############
# Plot MultiHiCCompare Results
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
