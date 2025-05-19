library(glue)
library(tidyverse)
library(furrr)
library(ggplot2)
library(ggpubr)
##############
# Load contacts + regions annotated
load_annotated_contacts <- function(
    regions.of.interest,
    ...){
    # Load all contacts
    load_mcool_files(
        pattern='*.NSC.*.mapq_30.1000.mcool',
        range1s='chr16',
        resolutions=c(100000)
    ) %>%
    mutate(is.Merged=grepl('Merged', Sample.ID)) %>%
    separate_wider_delim(
        Sample.ID,
        delim=fixed('.'),
        cols_remove=FALSE,
        names=c(
            'Edit',
            'Genotype',
            'SampleNumber',
            'Celltype'
        )
    ) %>% 
    annotate_contact_regions(
        regions.of.interest=regions.of.interest,
        most_specific_only=TRUE
    )
}
##############
# Plotting 
plot_contacts_regions_boxplot <- function(
    plot.df,
    bin_col='region',
    count_col='IF',
    color_col='region.color',
    facet_row='ReplicateNum',
    facet_col='Genotype',
    test.method='t.test',
    p.adjust.method='none',
    independent='none',
    scales='free_y',
    p.max=1,
    tip.length=0.05,
    n_pos=0,
    n_size=2,
    p_size=2,
    test_offset_y=3,
    expansion=c(0.01, 0.01, 0.1, 0.01),
    ...){
    stat.df <- 
        plot.df %>% 
        group_by(across(all_of(c(facet_row, facet_col)))) %>% 
        nest() %>% 
        rowwise() %>% 
        # Compute stats comparing means between bin_col (x-axis) values
        mutate(
            stats=
                compare_means(
                    formula(glue('{count_col} ~ {bin_col}')),
                    method=test.method,
                    p.adjust.method=p.adjust.method,
                    data=data
                ) %>%
                list()
        ) %>% 
        select(-c(data)) %>%
        unnest(stats) %>% 
        filter(p < p.max) %>% 
        # Make more legible label
        mutate(
            p.label=
                ifelse(
                    p.adjust.method == 'none',
                    glue('p={formatC(p, format="e", digits=2)}'),
                    glue('{p.adjust.method} p={formatC(p.adj, format="e", digits=2)}')
                )
        ) %>% 
        # Set positions + spacing of p-value text
        left_join(
            plot.df %>% 
            group_by(across(c(facet_row, facet_col))) %>% 
            summarize(
                y.position=
                    min(!!sym(count_col)) + 0.95 * (max(!!sym(count_col)) - min(!!sym(count_col)))
            ),
            by=c(facet_row, facet_col)
        ) %>%
        group_by(across(c(facet_row, facet_col))) %>% 
        mutate(y.position=max(y.position) + test_offset_y * (row_number())) %>% 
        ungroup()
    plot.df %>% 
    ggplot(
        aes(
            x=.data[[bin_col]],
            y=.data[[count_col]],
        )
    ) +
    geom_boxplot(
        aes(color=.data[[color_col]]),
        outlier.size=0.5
    ) +
    scale_color_identity() +
    stat_pvalue_manual(
        data=stat.df,
        label='p.label',
        size=p_size,
        tip.length=tip.length
    ) +
    geom_text(
        data=
            plot.df %>% 
            group_by(across(c(bin_col, facet_row, facet_col))) %>%
            summarize(
                y.position=n_pos,
                n=glue('n={n()}')
            ),
        aes(
            y=y.position,
            label=n
        ),
        size=n_size,
        hjust=0.5,
        vjust=0.5
    ) +
    geom_hline(yintercept=0, linetype='solid', color='black', linewidth=0.5) +
    scale_y_continuous(expand=expansion) +
    labs(y='HiC-Contacts') +
    facet_grid2(
        rows=vars(!!sym(facet_row)), 
        cols=vars(!!sym(facet_col)),
        scales=scales,
        independent=independent
    ) +
    theme(
        legend.position='none',
        axis.text.x=element_text(angle=25, hjust=1),
        axis.title.x=element_blank(),
        ...
    ) +
    add_ggtheme()
}
