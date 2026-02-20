###################################################
# Dependencies
###################################################
# library(glue)
# library(tidyverse)
# library(furrr)
# library(ggplot2)
# library(ggpubr)
library(hictkR)

###################################################
# Region comparisonRegion comparison
###################################################
list_all_region_comparisons <- function(
    huper.params.df,
    region.comparisons){
    # List region information for all regions
    regions.of.interest <- 
        GENOMIC_REGIONS %>%
        select(-c(region.group)) %>% 
        rename_with(~ str_remove(.x, '^region.')) %>%
        select(
            -c(
                chr,
                start,
                end
            )
        )
    # Get info for all regions to compare + background effect
    region.comaprisons.df <- 
        region.comparisons %>% 
        left_join(
            regions.of.interest,
            relationship='many-to-many',
            by=join_by(region.P1 == region)
        ) %>% 
        left_join(
            regions.of.interest,
            relationship='many-to-many',
            suffix=c('.P1', '.P2'),
            by=join_by(region.P2 == region)
        ) %>% 
        left_join(
            regions.of.interest %>%
            rename_with(~ str_replace(.x, '$', '.background')),
            relationship='many-to-many',
            by=join_by(region.background)
        )
    # List all hic matrices to compare regions within
    samples.df <- 
        list_mcool_files(only_use_included_samples=FALSE) %>%
        # filter(isMerged) %>%
        select(filepath, SampleID, isMerged)
    # Join everything together
    samples.df %>% 
        cross_join(hyper.params.df) %>% 
        cross_join(region.comaprisons.df)
}

load_contacts_to_compare <- function(
    region.P1, UCSC.P1,
    region.P2, UCSC.P2,
    UCSC.background,
    SampleID, filepath,
    resolution,
    normalization,
    ...){
    # Get contacts from the pair of regions of interest (RoI) being tested
    regions.contacts <-
        load_mcool_file(
            filepath,
            resolution=resolution,
            normalization=normalization,
            range1=UCSC.P1,
            range2=UCSC.P2
        ) %>%
        mutate(dist=range2 - range1)
    # Get distance range of RoI to filter background contacts
    dist.range <- 
        regions.contacts %>%
        summarize(
            max=max(dist),
            min=min(dist)
        )
    # Get background contact data to compare against RoI
    background.contacts <- 
        load_mcool_file(
            filepath,
            resolution=resolution,
            normalization=normalization,
            range1=UCSC.background
        ) %>%
        # Only keep distance-matched controls
        mutate(dist=range2 - range1) %>%
        filter(
            dist >= dist.range$min[[1]],
            dist <= dist.range$max[[1]]
        ) %>%
        # Remove regions of interest from background
        anti_join(
            regions.contacts %>%
            select(range1, range2),
            by=join_by(range1, range2)
        )
    # return al lcontact data to compare
    bind_rows(
        regions.contacts %>% add_column(group='region.of.interest'),
        background.contacts %>% add_column(group='distance.matched.background')
    )
}

###################################################
# Plotting 
###################################################
plot_contacts_regions_boxplot <- function(
    plot.df,
    resolution,
    normalization,
    test.method,
    results_file,
    stat.var,
    y.var,
    facet.row,
    facet.col,
    label.x.npc="left",
    label.y.npc="top",
    width=8,
    height=8,
    scales='free_y',
    n_size=4,
    ylim=TRUE,
    exclude_background=FALSE,
    ...){
    paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', collapse='; ')
    # row_index=1; isMerged=tmp$isMerged[[row_index]]; test.method=tmp$test.method[[row_index]]; resolution=tmp$resolution[[row_index]]; normalization=tmp$normalization[[row_index]]; stat.var=tmp$stat.var[[row_index]]; facet.row=tmp$facet.row[[row_index]]; comparison.group=tmp$comparison.group[[row_index]]; resolution.label=tmp$resolution.label[[row_index]]; region.comparison=tmp$region.comparison[[row_index]]; results_file=tmp$results_file[[row_index]]; plot.df=tmp$plot.df[[row_index]]
    # label.x.npc="left"; label.y.npc="top"; width=8; height=8; scales='free_y'; n_size=4
    # Load contacts to compare
    data.df <- 
        plot.df %>%
        mutate(
            contacts.df=
                pmap(
                    .l=.,
                    # .f=compare_interaction_vs_background,
                    .f=load_contacts_to_compare,
                    resolution=resolution,
                    normalization=normalization,
                    test.method=test.method,
                    .progress=FALSE
                )
        ) %>% 
        unnest(contacts.df) %>% 
        separate_wider_delim(
            SampleID,
            delim='.',
            names=c(NA, 'Celltype', 'Genotype', NA, NA),
            cols_remove=FALSE
        ) %>% 
        {
            if (exclude_background) {
                filter(., group != 'distance.matched.background')
            } else {
                .
            }
        } %>% 
        mutate(
            group=
                case_when(
                    group == 'region.of.interest'          ~ glue('{region.P2} - {region.P1}'),
                    group == 'distance.matched.background' ~ glue('[{region.background}] Distance Matched Control'),
                    TRUE                                   ~ NA
                )
        )
    # List stats comparisons to test
    my_comparisons <- 
        data.df[[stat.var]] %>%
        unique() %>%
        combn(m=2, simplify=FALSE)
    # Make boxplot of counts
    plot.obj <- 
        ggplot(
            data=data.df,
            aes(
                x=.data[[stat.var]],
                y=.data[[y.var]],
                fill=.data[[stat.var]]
                # x=group,
                # y=IF,
                # fill=group
            )
        ) + 
        geom_boxplot(outlier.size=0.4) +
        # geom_boxplot(outlier.shape=NA) +
        facet_grid(
            rows=vars(!!sym(facet.row)),
            cols=vars(!!sym(facet.col)),
            # rows=vars(Genotype),
            # cols=vars(Celltype),
            scales='free_y'
        ) +
        # compute t.test between groups
        stat_compare_means(
            comparisons=my_comparisons,
            size=n_size,
            label.x.npc=label.x.npc,
            label.y.npc=label.y.npc,
            method=test.method,
            data=data.df
        ) +
        # list number of bins being compared
        geom_text(
            data=
                data.df %>% 
                filter(IF > 0) %>% 
                # group_by(SampleID, Celltype, Genotype, group) %>%
                group_by(across(c('SampleID', facet.row, facet.col, stat.var))) %>%
                summarize(
                    n=dplyr::n(),
                    # y.position=150
                    y.position=quantile(IF, 0.82)
                ) %>% 
                mutate(label=glue('# Bin Pairs={n}')),
            aes(
                y=y.position,
                color=.data[[stat.var]],
                label=label
            ),
            show.legend=FALSE,
            size=n_size,
            hjust=0.5,
            vjust=0.5
        ) +
        # scale_y_continuous(expand=expansion) +
        labs(y='HiC-Contacts') +
        theme(
            legend.position='none',
            axis.text.x=element_text(angle=35, hjust=1),
            axis.title.x=element_blank(),
            # ...
        ) +
        make_ggtheme()
    if (ylim){
        plot.obj <- plot.obj + coord_cartesian(ylim=quantile(data.df[[y.var]], c(0.0, 0.98)))
    }
    ggsave(
        filename=results_file,
        plot=plot.obj,
        width=width,
        height=height,
        units='in'
    )
}

DEP_plot_contacts_regions_boxplot <- function(
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
        # group_by(across(all_of(c(facet_row, facet_col)))) %>% 
        # nest() %>% 
        # rowwise() %>% 
        # # Compute stats comparing means between bin_col (x-axis) values
        # mutate(
        #     stats=
        #         compare_means(
        #             formula(glue('{count_col} ~ {bin_col}')),
        #             method=test.method,
        #             p.adjust.method=p.adjust.method,
        #             data=data
        #         ) %>%
        #         list()
        # ) %>% 
        # select(-c(data)) %>%
        # unnest(stats) %>% 
        # filter(p < p.max) %>% 
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
    # geom_hline(yintercept=0, linetype='solid', color='black', linewidth=0.1) +
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
    scale_y_continuous(
        expand=expansion,
        breaks=c(0, 25, 50, 75)
    ) +
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

