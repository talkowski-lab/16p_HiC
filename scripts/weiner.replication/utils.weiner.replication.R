###################################################
# Dependencies
###################################################
library(glue)
library(tidyverse)
library(furrr)
library(ggplot2)
library(ggpubr)

###################################################
# Annotate contacts with specified regions
###################################################
annotate_contact_region_pairs <- function(
    contacts.df,
    regions.of.interest,
    most_specific_only=TRUE,
    ...){
    # regions of interest is a tibble where each row is a specific contiguous genomic range
    # see GENOMIC_REGIONS variable for example
    # For each bin (1, left and 2, right) in each pair of bins annotate the region is belongs to
    region.pairs.of.interest <- 
        regions.of.interest %>%
        # all possible pairs of regions to annotate
        # each bin-pair will be assigned 1 region pair, even if regions overlap/contain eachother
        join_all_rows(
            .,
            {.},
            suffix=c('1','2'),
            relationship='many-to-many',
        ) %>% 
        filter(region.chr1 == region.chr2) %>% 
        # factor order so most specific region interactions have the lowest factor level
        mutate(
            region=
                pmap_chr(
                    list(region1, region2),
                    ~ paste(sort(c(...)), collapse=' vs ')
                ),
            region.dist=region.dist1 + region.dist2,
            region=fct_reorder(region, region.dist) 
        )
    # contacts.df should be a tibble where each row is a pair of genomic bins i.e. output of
    # the load_mcool_file(s) functions
    contacts.df %>% 
    # Now annotate all region-pairs that each bin-pair belongs to 
    # (bin1 in region1 and bin2 in region2)
    full_join(
        region.pairs.of.interest,
        relationship='many-to-many',
        by=
            join_by(
                chr == region.chr1,
                between(
                    range1,
                    region.start1,
                    region.end1
                ),
                between(
                    range2,
                    region.start2,
                    region.end2
                )
            )
    ) %>% 
    # Pick most specific (smallest) annotated region for each bin if specified
    {
        if (most_specific_only) {
            group_by(
                .,
                SampleID,
                chr,
                range1,
                range2
            ) %>% 
            slice_min(
                region.dist,
                n=1,
                with_ties=FALSE
            ) %>% 
            ungroup()
        } else {
            .
        }
    } %>% 
    # Clean up 
    select(-c(starts_with('region.')))
}

format_annotated_contacts <- function(contacts.df) {
    contacts.df %>% 
    # Filter for relevant contacts
    mutate(contact.dist=range2 - range1) %>% 
    filter(
        contact.dist > dist.min,
        contact.dist < dist.max
    ) %>% 
    # Number samples for faceting plots
    nest(
        data=
            -c(
                Edit,
                Celltype,
                Genotype,
                is.Merged, 
                Sample.ID
            )
    ) %>% 
    arrange(
        is.Merged,
        Genotype
    ) %>% 
    # group_by(
    #     Edit,
    #     Celltype,
    #     Genotype
    # ) %>% 
    # mutate(ReplicateNum=as.character(row_number())) %>%
    # ungroup() %>% 
    mutate(
        ReplicateNum=
            ifelse(
                grepl('Merged', Sample.ID),
                'Merged',
                glue('Bio Replicate {ReplicateNum}')
            )
    ) %>% 
    unnest(data) %>% 
    # Plot visual details
    left_join(
        colors.and.titles,
        by='region'
    ) %>% 
    mutate(
        region.color=
            ifelse(
                is.na(region.color),
                '#808080', 
                region.color
            ),
        Genotype=
            factor(
                Genotype,
                levels=
                    c(
                        'WT',
                        'DEL',
                        'DUP'
                    )
            )
    ) %>% 
    # Now produce results with different minimum IF thresholds
    lapply(
        0:5,
        function(threshold, df) {
            df %>% 
            filter(IF > threshold) %>% 
            add_column(IF.Threshold=glue('IF > {threshold}'))
        },
        df=. 
    ) %>%
    bind_rows()
}

###################################################
# Plotting 
###################################################
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

