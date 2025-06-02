library(tidyverse)
library(stringi)
# library(HiContacts)
library(glue)
library(furrr)
library(ggplot2)
library(viridis)
library(ggh4x)
# library(ggtext)
library(cowplot)
###############
# Load TAD annotation data
load_hiTAD_DIs <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'start',
                'end',
                'DI'
            )
    )
}

load_arrowhead_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        skip=2,
        show_col_types=FALSE,
        col_select=c(1,2,3,12,13,14),
        # col_select=c(1,2,3,4,5,6,12,13,14),
        col_names=
            c(
                'chr',
                'start',
                'end',
                'B.chr',
                'B.start',
                'B.end',
                'name',
                'score',
                'strand1',
                'strand2',
                'color',
                'corner.score',
                'uVarScore',
                'lVarScore',
                'upSign',
                'loSign'
            )
    )
}

load_hiTAD_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'TAD.start',
                'TAD.end'
            )
    )
}

load_cooltools_TAD_annotation <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        progress=FALSE,
        show_col_types=FALSE
    ) %>%
    filter(!is_bad_bin) %>% 
    filter(if_any(starts_with('is_boundary'), ~ .x)) %>% 
    # select(starts_with('is_boundary')) %>% group_by(across(everything())) %>% count()
    rename(
        'chr'=chrom,
        'TAD.start'=start,
        'TAD.end'=end
    ) %>% 
    pivot_longer(
        ends_with('00'),
        names_to='TAD.stat',
        values_to='value'
    ) %>%
    select(-c(region, is_bad_bin)) %>% 
    mutate(TAD.stat=stri_reverse(TAD.stat)) %>% 
    separate_wider_delim(
        TAD.stat,
        delim='_',
        names=c(
            'window.size',
            'stat'
        ),
        too_many='merge'
    ) %>%
    mutate(
        across(
            c(stat, window.size),
            stri_reverse
        ),
        window.size=as.integer(window.size)
    ) %>% 
    filter(stat %in% c('boundary_strength', 'is_boundary', 'n_valid_pixels')) %>%
    pivot_wider(
        names_from=stat,
        values_from=value
    ) %>%
    mutate(is_boundary=as.logical(is_boundary))
}

load_TAD_annotation <- function(
    filepath,
    method,
    ...){
    if (method == 'hiTAD') {
        load_hiTAD_TAD_annotation(filepath)
    } else if (method == 'hiTAD-DIs') {
        load_hiTAD_DIs(filepath)
    } else if (method == 'cooltools') {
        load_cooltools_TAD_annotation(filepath)
    } else if (method == 'arrowhead') {
        load_arrowhead_TAD_annotation(filepath)
    } else {
        message(glue('Invalid method {method}'))
    }
}

load_all_TAD_annotations <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-TAD.tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    process_matrix_name() %>% 
    mutate(
        TADs=
            pmap(
                .,
                load_TAD_annotation,
                .progress=TRUE
            )
    ) %>%
    select(-c(filepath)) #%>% unnest(TADs)
}
###############
# Compute stuff

###############
# Plot stuff
plot_nTADs_heatmap <- function(
    plot.df,
    x_var='',
    y_var='',
    fill_var='',
    facet_row=NULL,
    facet_col=NULL,
    scales='fixed',
    ...){
    figure <- 
        plot.df %>% 
        ggplot(
            aes(
                x=.data[[x_var]],
                y=.data[[y_var]],
                fill=.data[[fill_var]],
            )
        ) +
        geom_tile() +
        scale_fill_viridis(discrete=FALSE) +
        # labs(title=sprintf('# of TADs TADs @ %sKb', .x$Resolution / 1000)) +
        theme(
            legend.position='top',
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1)
        ) +
        add_ggtheme()
    figure <- 
        add_faceting(
            figure,
            facet_col=facet_col,
            facet_row=facet_row,
            scales=scales
        )
    figure
}

plot_TADlengths_boxplot <- function(
    plot.df,
    x_var='',
    y_var='',
    color_var='',
    facet_row=NULL,
    facet_col=NULL,
    scales='fixed',
    ...){
    figure <- 
        plot.df %>% 
        ggplot(
            aes(
                y=.data[[y_var]],
                x=.data[[x_bar]],
                color=.data[[color_var]],
            )
        ) +
        geom_boxplot() +
        theme(
            axis.title.x=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_markdown()
        ) +
        add_ggtheme()
    figure <- 
        add_faceting(
            figure,
            facet_col=facet_row,
            facet_row=facet_col,
            scales=scales
        )
    figure
}
