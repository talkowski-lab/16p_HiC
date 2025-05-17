library(glue)
library(tidyverse)
library(furrr)
library(ggplot2)
library(ggpubr)
##############
# Plotting utils
plot_elise_boxplot <- function(
    plot.df,
    p_size=NA,
    n_size=NA,
    n_pos=0,
    independent='none',
    scales='free_y',
    test_offset_y=2,
    expansion=c(0.01, 0.01, 0.1, 0.01),
    ...){
    group.colors <- 
        c(
            '16all.Distance.Matched.Control'='#43090c',
            '16p.Distance.Matched.Control'='#e32528',
            '16q.Distance.Matched.Control'= '#871218',
            '16p11.2 CNV - Telomere'='#526ab1'
        )
    stat.df <- 
        plot.df %>% 
        group_split(Replicate.ID, Genotype) %>% 
        lapply(
            function(df) {
                compare_means(
                    formula(glue('count ~ region')),
                    method='t.test',
                    p.adjust.method='BH',
                    data=df
                ) %>% 
                group_by(group1, group2) %>% 
                add_column(
                    Genotype=df$Genotype[[1]],
                    Replicate.ID=df$Replicate.ID[[1]]
                )
            }
        ) %>%
        bind_rows() %>% 
        left_join(
            plot.df %>% 
            group_by(Replicate.ID, Genotype) %>% 
            summarize(y.position=min(count) + 0.75 * (max(count) - min(count))),
            by=
                join_by(
                    Genotype,
                    Replicate.ID,
                )
        ) %>%
        group_by(Replicate.ID, Genotype) %>% 
        mutate(y.position=max(y.position) + test_offset_y * (row_number())) %>% 
        ungroup()
    plot.df %>% 
        ggplot(
            aes(
                x=region,
                y=count,
                color=region
            )
        ) +
        geom_boxplot(outlier.size=0.5) +
        scale_y_continuous(
            expand=expansion
        ) +
        stat_pvalue_manual(
            stat.df,
            label='T-test p={formatC(p.adj, format="e", digits=2)}',
            size=p_size,
            tip.length=0
        ) +
        geom_text(
            data=
                plot.df %>% 
                group_by(Genotype, Replicate.ID, region) %>%
                summarize(
                    # y.position=max(.data[[y_val]]),
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
        scale_color_manual(values=group.colors) +
        labs(y='HiC-Contacts') +
        facet_grid2(
            rows=vars(Replicate.ID),
            cols=vars(Genotype),
            scales=scales,
            independent=independent
        ) +
        theme(
            legend.position='none',
            axis.text.x=element_text(angle=25, hjust=1),
            axis.title.x=element_blank()
        ) +
        add_ggtheme()
}

plot_wrapper <- function(
    interactions,
    min_count,
    plot_file,
    independent='y',
    scales='free_y',
    p_size=2.5,
    n_size=2.75,
    test_offset_y=5,
    expansion=1e-2 * c(4, 1, 10, 1),
    height=10,
    width=10,
    ...){
    interactions %>% 
    filter(count > min_count) %>% 
    plot_elise_boxplot(
        independent=independent,
        scales=scales,
        p_size=p_size,
        n_size=n_size,
        test_offset_y=test_offset_y,
        expansion=expansion
    )
}
