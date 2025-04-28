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

make_tab_per_group <- function(
    plot.df, 
    group_col,
    plot_fnc,
    col_header=TRUE,
    max_lvl=3,
    ...){

    cat("\n\n")
    max_header <- rep('#', max_lvl) %>% paste0(collapse="")
    cat(max_header, ifelse(col_header, group_col, ""), "{.tabset}", "\n\n")
    tab_header <- paste0(max_header, '#')
    group_values <- 
        plot.df[[group_col]] %>% 
        as.factor() %>% 
        droplevels() %>%
        levels()
    for (group_value in group_values){
        cat(tab_header, group_value, "\n\n")
        figure <- 
            plot.df %>%
            filter(get({{group_col}}) == group_value) %>%
            plot_fnc(...)
        print(figure)
        cat("\n\n")
    }
}
