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

make_pair_qc_df <- function(
    pair.stats.list,
    cols_to_keep=c(),
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
        across(
            c('total.unique.contacts', cols_to_keep),
            ~ unique(.x)
        )
    ) %>% 
    ungroup() %>% 
    mutate(frac=(value / total.unique.contacts))
}

make_summary_stats_table <- function(
    min.resolutions,
    pair.stats.list,
    pairs.plot.df,
    ...){
    # Start with coverage data
    min.resolutions %>% 
    select(
        SampleID,
        ReadFilter,
        # metric,
        # chr,
        `Min. Viable Resolution`,
        `# Total Bins`,
        `# Bins > 0`,
        `% Bins > 0`,
        `# Bins >= 1K`,
        `% Bins >= 1K`
    ) %>% 
    rename('Min. MAPQ Threshold'=ReadFilter) %>% 
    mutate(across(where(is.numeric), as.character)) %>% 
    pivot_longer(
        -c(SampleID),
        names_to='metric',
        values_to='value'
    ) %>%
    # add pair categories
    bind_rows(
        pairs.plot.df %>%
        select(
            SampleID,
            Category,
            value
        ) %>%
        rename(
            'metric'=Category,
            'value'=value
        ) %>%
        mutate(
            metric=paste0('pairs.', metric),
            value=as.character(value)
        )
    ) %>% 
        # distinct(SampleID) %>% print(n=Inf)
        # {.} %>% head(2) %>% t()
    # Add stats per matrix
    bind_rows(
        pair.stats.list$general.stats %>%
        select(
            SampleID,
            stat,
            value
        ) %>%
        mutate(value=as.character(value)) %>% 
        rename('metric'=stat)
    ) %>% 
    pivot_wider(
        names_from=metric,
        values_from=value
    ) %>%
    arrange(SampleID) %>%
        pivot_longer(-c(SampleID)) %>%
        group_by(SampleID) %>% 
        summarize(n.missing=sum(is.na(value))) %>% 
        arrange(desc(n.missing))
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
plot_barplot_with_threshold_lines <- function(
    plot.df,
    lines.df,
    y_lower=0.15,
    y_upper=0.40,
    line_bleed=0.65,
    linewidth=0.4,
    linetype='dashed',
    ...){
    # lines.df <- 
    #     tribble(
    #         ~Category, ~y_lower, ~y_upper,
    #         'long range', 0.15, 0.40
    #     )
    lines.variable <- 
        lines.df %>% select(-c(y_lower, y_upper)) %>% colnames() %>% unlist() 
    # Define where to draw horizontal threshold lines and how long across bars
    lines.df %>% 
    mutate(
        !!lines.variable := 
                factor(
                    !!sym(lines.variable),
                    levels=levels(plot.df[[lines.variable]])
                ),
        x_line=as.numeric(!!sym(lines.variable)) - line_bleed,
        xend_line=as.numeric(!!sym(lines.variable)) + line_bleed
    ) %>% 
    # add this to ploting data
    right_join(
        plot.df,
        by=lines.variable
    ) %>% 
        # {.}
    # Draw regular barplot 
    plot_barplot(...) +
    # Draw threshold lines on top of the barplot with parameters set above
    geom_segment(
        aes(
            x=x_line,
            xend=xend_line
        ),
        y=y_lower,
        yend=y_lower,
        color='red',
        linewidth=linewidth,
        linetype=linetype
    ) +
    geom_segment(
        aes(
            x=x_line,
            xend=xend_line
        ),
        y=y_upper,
        yend=y_upper,
        color='green',
        linewidth=linewidth,
        linetype=linetype
    )
}

plot_pair.orientation_lineplot <- function(
    plot.df,
    facet_row,
    legend.position.inside=c(0.95, 0.25), # c(0,0) bottom left, c(1,1) top-right.
    ...){
    plot.df %>% 
    ggplot(
        aes(
            x=range.start,
            y=value,
        )
    ) + 
    # geom_line(aes(linetype=orientation, color=SampleID)) +
    geom_line(aes(color=orientation)) +
    # geom_text(aes(label=SampleID), x=1.5, y=1.8) +
    scale_x_log10(
        guide='axis_logticks',
        labels=
            label_log(
                base=10,
                signed=FALSE
            )
    ) +
    # scale_x_continuous(
    #     labels=
    #         label_bytes(
    #             units="auto_si",
    #             accuracy=0.1
    #         )
    # ) +
    scale_y_log10(
        guide='axis_logticks',
        labels=
            label_log(
                base=10,
                signed=FALSE
            )
    ) +
    # facet_wrap( ~ SampleID, ncol=3) +
    facet_grid(
        rows=vars(!!sym(facet_row)),
        scales='fixed'
    ) +
    labs(
        x='Pair Distance',
        y='Raw Contacts'
    ) +
    make_ggtheme(
        legend.position='inside',
        legend.position.inside=legend.position.inside,
        ...
    )
}

plot_cistrans_chromosome_heatmap <- function(
    plot.df,
    fill_var,
    facet_group,
    scales='fixed',
    legend.position='top',
    axis.text.x=element_text(angle=45, hjust=1),
    ...){
    plot.df %>%
    plot_heatmap(
        x_var='chr1',
        y_var='chr2',
        fill_var=fill_var, 
        facet_row=NULL,
        facet_col=NULL,
        facet_group=facet_group,
        legend.position=legend.position,
        axis.text.x=axis.text.x,
        scales=scales,
        ...
    ) +
    scale_fill_gradient(
        low='grey90',
        high='red'
    ) +
    theme(axis.title=element_blank())
}
