library(glue)
library(tidyverse)
library(furrr)
library(ggplot2)
library(ggpubr)
library(hictkR)
###############
# Dirs
BASE_DIR="/data/talkowski/Samples/16p_HiC"
BASE_DIR="/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p"
source(glue('{BASE_DIR}/scripts/utils.data.R'))
RESULTS_DIR=glue('{ELISE_RECREATION_DIR}/results')
PLOT_DIR=glue('{ELISE_RECREATION_DIR}/plots')
###############
# Params
RESOLUTIONS=c(100000)
NORMALIZATIONS=c('cis.ice', 'cis.vc', 'NONE')
# Important Regions
CHR16_ELISE_REGION='chr16:29,488,679-30,188,679'
CHR16_TELOMERE_REGION='chr16:0-5,149,999'
CHR16_CONTROL_REGIONS=
    list(
        '16p11.2'='chr16:24,288,679-30,188,679',
        '16p'=    'chr16:0-36,800,000',
        '16q'=    'chr16:36,800,001-90,338,345',
        '16.all'= 'chr16:0-90,338,345'
    )
###############
# Fetch contacts region 
fetch_contacts <- function(
    mcool_file,
    resolution,
    normalization,
    ...){
    # hic.obj=hic_matrices$matrix[[1]]; SampleID=hic_matrices$SampleID[[1]]
    interactions.file <- mcool_file %>% File(resolution=resolution)
    cnv.interactions <- 
        interactions.file %>%
        fetch(
            range1=CHR16_TELOMERE_REGION,
            range2=CHR16_ELISE_REGION,
            normalization=normalization,
            join=TRUE,
            query_type='UCSC',
            type='df'
        ) %>% 
        as_tibble() %>% 
        add_column(region='16p11.2 CNV - Telomere')
    dist.range <- 
        cnv.interactions %>%
        mutate(dist=start2 - start1) %>%
        {list('min'=min(.$dist), 'max'=max(.$dist))}
    telomere <- 
        CHR16_TELOMERE_REGION %>% 
        str_remove('^chr16:') %>% 
        str_remove_all(',') %>% 
        str_split('-') %>% 
        first() 
    # Get Telomere+CNV boundaries to not double-count them in controls
    telomere.start <- as.integer(telomere[[1]])
    telomere.end <- as.integer(telomere[[2]])
    cnv <- 
        CHR16_ELISE_REGION %>% 
        str_remove('chr16:') %>% 
        str_remove_all(',') %>% 
        str_split('-') %>% 
        first()
    cnv.start <- as.integer(cnv[[1]])
    cnv.end <- as.integer(cnv[[2]])
    # Get sets of control interactions to compare against i.e.
    # all contacts within and distance range as Telomere-CNV contacts
    CHR16_CONTROL_REGIONS %>% 
    as_tibble() %>%
    # First get contacts within each region of interest
    pivot_longer(
        everything(),
        names_to='region',
        values_to='coords'
    ) %>%
    # Now find all contacts within specified distance range
    mutate(
        region=glue('{region}.Distance.Matched.Control'),
        interactions=
            pmap(
                .l=.,
                .f=function(interactions.file, coords, ...){
                    interactions.file %>%
                    fetch(
                        range1=coords,
                        normalization=normalization,
                        join=TRUE,
                        query_type='UCSC',
                        type='df'
                    )
                },
                interactions.file=interactions.file
            )
    ) %>% 
    unnest(interactions) %>% 
    mutate(dist=start2 - start1) %>%
    dplyr::filter(between(dist, dist.range$min, dist.range$max)) %>%
    # Remove any contacts spanning regions of interest
    filter(
        !(
            between(start1, telomere.start, telomere.end) &
            between(start2, cnv.start, cnv.end)
        )
    ) %>% 
    select(-dist) %>% 
    # add cnv contacts back and return 
    bind_rows(cnv.interactions) %>%
    mutate(
        region=
            factor(
                region,
                levels=
                    CHR16_CONTROL_REGIONS %>% 
                    names() %>% 
                    paste0('.Distance.Matched.Control') %>% 
                    sort() %>% 
                    c('16p11.2 CNV - Telomere', .)
            )
    ) %>% 
    select(-c(coords, end1, end2))
}
# Import contact data
interactions.df <- 
    load_16p_sample_metadata() %>% 
    expand_grid(
        normalization=NORMALIZATIONS,
        resolution=RESOLUTIONS,
    ) %>%
    mutate(
        interactions=
            pmap(
                .l=.,
                .f=fetch_contacts,
                .progress=TRUE
            )
    )
###############
# Make boxplots
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
    ) %>%
    ggsave(
        plot_file,
        .,
        width=width,
        height=height,
        units='in'
    )
}
plot.df <- 
    interactions.df %>% 
    left_join(
        expand_grid(
            normalization='NONE',
            min_count=1:5
        ),
        by='normalization',
        relationship='many-to-many'
    ) %>% 
    mutate(min_count=ifelse(is.na(min_count), -Inf, min_count)) %>% 
    unnest(interactions) %>% 
    # mutate(count=ifelse(is.nan(count), 0, count)) %>% 
    select(-c(chrom2, mcool_file)) %>% 
    mutate(plot_file=glue('{PLOT_DIR}/{normalization}/{min_count}/contact.boxplot.pdf'))
plot.df %>% select(-c(plot_file)) %>% head(2) %>% t()
plot.df %>%
    filter(resolution == 100000) %>%
    filter(chrom1 == 'chr16') %>%
    filter(normalization == 'cis.ice') %>% 
    group_by(Genotype, Replicate.ID) %>% 
    summarize(n=n(), total=sum(count, na.rm=TRUE), as_tibble_row(summary(count))) %>%
    ungroup() #%>% t()
# cnv.interactions
###############
# Plot IFs of ELISE vs CTRL regions 
plot.df %>%
    filter(region != '16.all.Distance.Matched.Control') %>%
    # filter(count < 200) %>% 
    nest(interactions=-c(normalization, min_count, plot_file)) %>%
    rowwise() %>%
    mutate(
        height=10,
        width=10,
        independent='y',
        scales='free_y',
        p_size=2.5,
        n_size=2.75,
        n_pos=ifelse(normalization == 'NONE', 0, 3e-4),
        test_offset_y=
            case_when(
                normalization == 'cis.ice' ~ 5e-4,
                normalization == 'NONE' ~ 2,
                normalization == 'cis.vc' ~ 10,
            ),
        expansion=
            ifelse(
                normalization == 'NONE',
                list(c(4e-2, 1e-2, 10e-2, 1e-2 )), 
                list(c(5e-3, 1e-4,  5e-2, 1e-4))
            )
    ) %>% 
    pmap(
        .l=.,
        .f=plot_wrapper,
    )
## Presentation version
plot.df %>% 
    filter(Replicate.ID == 'Merged') %>% 
    filter(normalization == 'NONE') %>% 
    filter(min_count == 1) %>% 
    filter(region %in% c('16p11.2 CNV - Telomere', '16p.Distance.Matched.Control')) %>%
    plot_elise_boxplot(
        independent='y',
        scales='free_y',
        y_val='count',
        p_size=2.5, n_size=2.75, test_offset_y=5, expansion=1e-2 * c(3, 1, 5, 1)
    ) %>%
    ggsave(
        glue('{PLOT_DIR}/presentation.boxplot.pdf'),
        .,
        width=10,
        height=6,
        units='in'
    )
