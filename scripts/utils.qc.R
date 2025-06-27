library(tidyverse)
library(glue)
library(purrr)
###############
# Load various QC data files/sets of files
load_pairtools_stats <- function(
    stats_file_suffix='.HiC.hg38.dedup.stats',
    ...){
    # Load all stats
    all.stats.df <- 
        PAIRS_DIR %>%
        list.files(
            full.names=FALSE,
            recursive=TRUE,
            pattern=glue('*{stats_file_suffix}$')
        ) %>%
        tibble(fileinfo=.) %>% 
        mutate(filepath=file.path(PAIRS_DIR, fileinfo)) %>% 
        mutate(Sample.ID=str_remove(basename(fileinfo), stats_file_suffix)) %>% 
        select(-c(fileinfo)) %>% 
        mutate(
            stats=
                purrr::pmap(
                    .l=.,
                    function(filepath, ...){
                        filepath %>% 
                        read_tsv(
                            show_col_types=FALSE,
                            col_names=c('stat', 'value')
                        ) %>%
                        filter(stat != 'summary/dist_freq_convergence/strands_w_max_convergence_dist') %>%
                        mutate(value=as.numeric(value))
                    }
                )
        ) %>%
        select(-c(filepath)) %>%
        unnest(stats)
    all.stats.df
    all.stats.df <- 
        all.stats.df %>% 
        separate_wider_delim(
            Sample.ID,
            delim=fixed('.'),
            names=c('Edit', 'Genotype', 'SampleNumber', 'Celltype'),
            cols_remove=FALSE
        ) %>% 
        group_by(Edit, Genotype, Celltype, stat) %>%
        summarize(value=sum(value)) %>%
        ungroup() %>% 
        mutate(Sample.ID=glue('{Edit}.{Genotype}.Merged.{Celltype}')) %>%
        select(-c(Edit, Genotype, Celltype)) %>% 
        bind_rows(all.stats.df)
    # Get total number of contacts per sample
    total.contacts.df <- 
        all.stats.df %>% 
        filter(!grepl('/', stat)) %>%
        filter(stat == 'total_nodups') %>% 
        select(Sample.ID, value) %>%
        dplyr::rename('total.unique.contacts'=value)
    # Matrix-wide stats
    general.stats.df <- 
        all.stats.df %>% 
        filter(!grepl('/', stat)) %>%
        mutate(
            category=
                case_when(
                    grepl('total', stat) ~ 'Reads', 
                    TRUE ~ 'Contacts'
                )
        ) %>% 
        left_join(
            total.contacts.df,
            by='Sample.ID'
        )
    # stats about contact distances
    dist.stats.df <- 
        all.stats.df %>% 
        filter(grepl('^dist_freq', stat)) %>%
        separate_wider_delim(
            stat,
            delim='/',
            names=
                c(
                    'interaction',
                    'range',
                    'orientation'
                )
        ) %>%
        mutate(range=str_remove(range, fixed('+'))) %>% 
        separate_wider_delim(
            range,
            delim='-',
            names=c('range.start', 'range.end'),
            too_few='align_start'
        ) %>%
        mutate(across(starts_with('range'), ~ as.integer(.x))) %>% 
        mutate(interaction='cis') %>%
        left_join(total.contacts.df, by='Sample.ID')
    # Contact Frequency Stats per Chromosome
    chr.stats.df <- 
        all.stats.df %>% 
        filter(grepl('^chrom_freq', stat)) %>%
        separate_wider_delim(
            stat,
            delim='/',
            names=
                c(
                    'stat',
                    'chr1',
                    'chr2'
                )
        ) %>% 
        left_join(total.contacts.df, by='Sample.ID')
    # Return named list of all stuff
    return(
        list(
            'general.stats'=general.stats.df,
            'distance.stats'=dist.stats.df,
            'chr.stats'=chr.stats.df
        )
    )
}

make_summary_stats_table <- function(
    coverage.summary,
    pair.stats.list,
    pair.plot.df,
    ...){
# Start with coverage data
    coverage.summary %>%
    filter(
        weight == 'raw',
        metric == 'cis',
        chr == 'genome.wide',
        ReadFilter == 'mapq_30',
        # bins.pct.covered >= 0.80
    ) %>% 
    # Pick the smallest resolution that is valid for each sample
    mutate(resolution.is.valid=bins.pct.covered >= 0.80) %>% 
    group_by(Sample.ID) %>% 
    slice_min(
        tibble(desc(resolution.is.valid), resolution),
        n=1,
        with_ties=FALSE
    ) %>% 
    ungroup() %>% 
    select(
        -c(
            resolution.is.valid,
            weight,
            metric,
            chr
        )
    ) %>%
    mutate(ReadFilter=as.integer(str_remove(ReadFilter, 'mapq_'))) %>% 
    pivot_longer(
        c(starts_with('bins'), resolution, ReadFilter),
        names_to='metric',
        values_to='n'
    ) %>%
    mutate(
        metric=
            case_when(
                metric == 'ReadFilter'       ~ 'Minimum MAPQ',
                metric == 'resolution'       ~ 'Minimum Usable Resolution',
                metric == 'bins.n.total'     ~ '# of Bins',
                metric == 'bins.n.nz'        ~ '# of Bins > 0 Contacts',
                metric == 'bins.pct.nz'      ~ '% Bins > 0 Contacts',
                metric == 'bins.n.covered'   ~ 'Bins >= 1K Contacts',
                metric == 'bins.pct.covered' ~ '% Bins >= 1K Contacts',
                TRUE ~ NA
            )
    ) %>% 
    filter(!is.na(metric)) %>% 
    # add pair categories
    bind_rows(
        pair.plot.df %>%
        select(
            Sample.ID,
            Category,
            value
        ) %>%
        rename(
            'metric'=Category,
            'n'=value
        ) %>%
        mutate(metric=paste0('pairs.', metric))
    ) %>% 
    # Add stats per matrix
    bind_rows(
        pair.stats.list$general.stats %>%
        select(
            Sample.ID,
            stat,
            value
        ) %>%
        rename(
            'metric'=stat,
            'n'=value
        )
    ) %>% 
    pivot_wider(
        names_from=metric,
        values_from=n
    ) %>%
    arrange(Sample.ID)
}

make_pair_qc_df <- function(
    pair.stats.list,
    ...){
    pair.stats.list$general.stats %>% 
    select(-c(category)) %>% 
    filter(stat %in% c('trans', 'cis')) %>% 
    dplyr::rename('interaction'=stat) %>% 
    add_column(range.start=0) %>% 
    bind_rows(pair.stats.list$distance.stats) %>% 
    mutate(
        range.end=ifelse(is.na(range.end), Inf, range.end),
        Category=
            case_when(
                orientation =='-+'   & range.end < 1001    ~ 'Self-Circle', 
                orientation =='+-'   & range.end < 1001    ~ 'Dangling-End', 
                range.start >= 20000                       ~ 'Long Range  >  20Kb',
                range.start >=  1000 & range.start < 20001 ~ 'Short Range <= 20Kb',
                range.end   <=  1000                       ~ 'Too Short   <=  1Kb',
                interaction == 'cis' & range.end == Inf    ~ 'Total Cis',
                interaction == 'trans'                     ~ 'Trans',
            )
    ) %>%
    mutate(
        Category=
            factor(
                Category,
                levels=
                    c(
                        'Self-Circle', 
                        'Dangling-End', 
                        'Too Short   <=  1Kb',
                        'Short Range <= 20Kb',
                        'Long Range  >  20Kb',
                        'Total Cis',
                        'Trans'
                    )
            )
    ) %>% 
    group_by(Sample.ID, Category) %>% 
    summarize(
        value=sum(value),
        # total per sample
        total.unique.contacts=unique(total.unique.contacts)
    ) %>% 
    ungroup() %>% 
    mutate(frac=(value / total.unique.contacts))
}

get_matrix_minimum_resolutions <- function(
    input_data,
    ...){
    {
        if (is_tibble(input_data)) {
            input_data
        } else {
            load_genome_coverage(input_data)
        }
    } %>% 
    group_by(across(-c(start, end, coverage))) %>% 
    summarize(
        bins.n.nz=sum(coverage > 0),
        bins.n.covered=sum(coverage > 1000),
        bins.n.total=n()
    ) %>%
    ungroup() %>% 
    # Calculate genomw-wide totals as well
    bind_rows(
        .,
        {.} %>% 
        group_by(
            across(
                -c(
                    chr,
                    starts_with('bins.n.')
                )
            )
        ) %>% 
        summarize(
            across(
                starts_with('bins.n.'),
                ~ sum(.x)
            )
        ) %>%
        ungroup() %>% 
        add_column(chr='genome.wide')
    ) %>%
    mutate(
        bins.pct.nz=bins.n.nz / bins.n.total,
        bins.pct.covered=bins.n.covered / bins.n.total
    )
}

get_all_matrix_minimum_resolutions <- function(
    resolutions=NULL,
    ...){
    # List input files (generated by cooltools coverage)
    parse_results_filelist(
        input_dir=COVERAGE_DIR,
        suffix='-coverage.tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    process_matrix_name() %>% 
    # Only need raw coverage
    filter(weight == 'raw') %>% 
    # Ignore results not matching param values
    { 
        if (is.null(resolutions)) { 
            . 
        } else { 
            filter(., resolution %in% resolutions) 
        } 
    } %>%
    mutate(resolution=as.integer(resolution)) %>%
    mutate(input_data=filepath) %>% 
    mutate(
        min.resolution.results=
            pmap(
                .l=.,
                .f=get_matrix_minimum_resolutions,
                .progress=TRUE
            )
    ) %>%
    unnest(min.resolution.results) %>%
    select(-c(input_data))
}
###############
# Plot Matrix QC stuff
plot_qc_barplot <- function(
    plot.df,
    y_val,
    scale_y_method,
    fill_col='Sample.ID',
    facet_col='Edit',
    facet_row='isMerged',
    expand=c(0.00, 0.00, 0.00, 0.00),
    scales='fixed',
    ...){
    figure <- 
        plot.df %>%
        ggplot(
            aes(
                x=Category,
                y=.data[[y_val]],
                fill=.data[[fill_col]]
            )
        ) +
        geom_col(position = "dodge") +
        guides(fill=guide_legend(ncol=1)) +
        facet_grid(
            rows=vars(!!sym(facet_row)),
            cols=vars(!!sym(facet_col)),
            scales=scales
        ) +
        labs(y=y_val) +
        theme(
            ...,
            axis.title.x=element_blank(),
            axis.text.x=
                element_text(
                    hjust=1,
                    angle=35
                )
        ) + 
        add_ggtheme()
    # scale axis
    figure %>% 
    scale_y_axis(
        mode=scale_y_method,
        expand=expand
    )
}

add_lines_qc_barplot <- function(
    plot.df,
    linewidth=0.4,
    linetype='dashed',
    ...){
    plot.df %>% 
    mutate(
        x_line=
            ifelse(
                grepl('long range', Category, ignore.case=TRUE), 
                as.numeric(Category) - 0.55,
                NA
            ),
        xend_line=
            ifelse(
                grepl('long range', Category, ignore.case=TRUE),
                as.numeric(Category) + 0.55,
                NA
            )
    ) %>% 
    plot_qc_barplot(...) +
    # Draw extra line on qc plot
    geom_segment(
        aes(
            x=x_line,
            xend=xend_line
        ),
        y=0.15, 
        yend=0.15,
        color='red',
        linewidth=linewidth,
        linetype=linetype
    ) +
    geom_segment(
        aes(
            x=x_line,
            xend=xend_line
        ),
        y=0.40,
        yend=0.40,
        color='green',
        linewidth=linewidth,
        linetype=linetype
    )
}

plot_pair.orientation_lineplot <- function(
    plot.df,
    facet_col,
    ...){
    plot.df %>% 
    ggplot(
        aes(
            x=range.start,
            y=value,
        )
    ) + 
    # geom_line(aes(linetype=orientation, color=Sample.ID)) +
    geom_line(aes(color=orientation)) +
    # geom_text(aes(label=Sample.ID), x=1.5, y=1.8) +
    scale_x_log10(labels=label_log()) +
    scale_y_log10(labels=label_log()) +
    # facet_wrap( ~ Sample.ID, ncol=3) +
    facet_grid(
        cols=vars(!!sym(facet_col)),
        scales='fixed'
    ) +
    labs(
        x='Pair Distance',
        y='Raw Contacts'
    ) +
    theme(
        legend.position='inside',
        legend.position.inside=c(0.95, 0.25), # c(0,0) bottom left, c(1,1) top-right.
        ...
    ) +
    add_ggtheme()
}

plot_cistrans_chromosome_heatmap <- function(
    plot.df,
    fill_var,
    facet_col,
    ...){
    plot.df %>%
    ggplot(
        aes(
            x=chr1,
            y=chr2
        )
    ) +
    geom_tile(
        aes(
            fill=.data[[fill_var]]
        )
    ) +
    scale_fill_gradient(
        low='grey90',
        high='red'
    ) +
    facet_grid(
        cols=vars(!!sym(facet_col)),
        scales='fixed'
    ) +
    theme(
        legend.position='top',
        axis.text.x=element_text(angle=45, hjust=1)
    ) +
    add_ggtheme()
}
