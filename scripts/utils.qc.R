###################################################
# Dependencies
###################################################
library(tidyverse)
library(glue)
library(purrr)

###################################################
# Load various QC data files/sets of files
###################################################
load_pairtools_stats <- function(
    stats_file_suffix='.hg38.dedup.stats',
    samples.to.include=NULL,
    SampleID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID'
        ),
    ...){
    # stats_file_suffix='.hg38.dedup.stats'; samples.to.include=NULL; SampleID.fields=c('Edit', 'Celltype', 'Genotype', 'CloneID', 'TechRepID', 'Batch')
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
        mutate(SampleID=str_remove(basename(fileinfo), stats_file_suffix)) %>% 
        {
            if (!is.null(samples.to.include)) {
                filter(., SampleID %in% samples.to.include)
            } else {
                .
            }
        } %>% 
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
                        # ignore stat that is not included in all outputs 
                        filter(stat != 'summary/dist_freq_convergence/strands_w_max_convergence_dist') %>%
                        mutate(value=as.numeric(value))
                    },
                    .progress=TRUE
                )
        ) %>%
        select(-c(filepath)) %>%
        unnest(stats)
    # Compute stats for merged matrices and add as new rows
    all.stats.df <- 
        all.stats.df %>% 
        get_info_from_SampleIDs(
            keep_id=FALSE,
            SampleID.fields=SampleID.fields
        ) %>% 
        group_by(across(all_of(c(setdiff(SampleID.fields, c('CloneID', 'TechRepID')), 'stat')))) %>%
        # group_by(Edit, Celltype, Genotype, stat) %>% 
        summarize(value=sum(value)) %>%
        ungroup() %>% 
        add_column(
            CloneID='Merged',
            TechRepID='Merged'
        ) %>% 
        mutate(SampleID=glue(paste0('{', SampleID.fields, '}', collapse='.'))) %>% 
        select(-all_of(SampleID.fields)) %>% 
        # mutate(SampleID=glue('{Edit}.{Celltype}.{Genotype}.{CloneID}.{TechRepID}')) %>% 
        # select(-c(Edit, Celltype, Genotype, CloneID, TechRepID)) %>% 
        bind_rows(all.stats.df)
    # all.stats.df %>% distinct(SampleID)
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
        mutate(interaction='cis')
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
        mutate(across(starts_with('chr'), ~ factor(.x, levels=CHROMOSOMES)))
    # Get total number of contacts per sample
    total.contacts.df <- 
        all.stats.df %>% 
        filter(!grepl('/', stat)) %>%
        filter(stat == 'total_nodups') %>% 
        select(SampleID, value) %>%
        dplyr::rename('total.unique.contacts'=value)
    # Return named list of all stuff + add total contacts in each sample for computing %
    return(
        list(
            'general.stats'=
                general.stats.df %>% 
                left_join(
                    total.contacts.df,
                    by='SampleID'
                ),
            'distance.stats'=
                dist.stats.df %>% 
                left_join(
                    total.contacts.df,
                    by='SampleID'
                ),
            'chr.stats'=
                chr.stats.df %>% 
                left_join(
                    total.contacts.df,
                    by='SampleID'
                )
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
    group_by(SampleID, Category) %>% 
    summarize(
        value=sum(value),
        across(
            all_of(c('total.unique.contacts', cols_to_keep)),
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

get_all_matrix_minimum_resolutions <- function(
    file.suffix='-coverage.tsv',
    resolutions=NULL,
    ...){
    # List input files (generated by cooltools coverage)
    parse_results_filelist(
        input_dir=COVERAGE_DIR,
        suffix=file.suffix,
        filename.column.name='MatrixID',
        param_delim='_',
    ) %>%
    mutate(MatrixID=str_remove(MatrixID, file.suffix)) %>% 
    get_info_from_MatrixIDs(
        ID_col='MatrixID',
        keep_id=FALSE,
        nest_col=NA,
        add_sample_id=TRUE,
    ) %>% 
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
    # load each coverage file and count % of bins > 0 and > 1000 contacts
    mutate(
        min.resolution.results=
            pmap(
                .l=.,
                .f=get_matrix_minimum_resolutions,
                .progress=TRUE
            )
    ) %>%
    unnest(min.resolution.results) %>%
    select(-c(filepath))
}

###################################################
# Genome Marginal Coverage
###################################################
list_all_coverage_data <- function(resolutions=NULL){
    # List input files (generated by cooltools coverage)
    COVERAGE_DIR %>% 
    parse_results_filelist(
        suffix='-coverage.tsv',
        filename.column.name='MatrixID',
        param_delim='_'
    ) %>%
    # mutate(MatrixID=str_remove(MatrixID, file.suffix)) %>% 
    get_info_from_MatrixIDs(
        suffix=file.suffix,
        keep_id=FALSE
    ) %>% 
    # only keep specified resoltuions
    { 
        if (is.null(resolutions)) { 
            . 
        } else { 
            filter(., resolution %in% resolutions) 
        } 
    }
}

load_all_coverage_data <- function(resolutions=NULL){
    list_all_coverage_data(resolutions=resolutions) %>% 
    # Only need raw cis coverage
    filter(weight == 'raw') %>% 
    filter(count.type == 'cis') %>% 
    # load each coverage file and count % of bins > 0 and > 1000 contacts
    mutate(
        coverage.data=
            pmap(
                .l=.,
                .f=load_genome_coverage,
                .progress=TRUE
            )
    ) %>%
    unnest(coverage.data) %>%
    select(-c(filepath))
}

compute_summary_stats <- function(df){
    df %>% 
    summarize(
        coverage.min=min(coverage),
        coverage.q25=quantile(coverage, 0.25),
        coverage.median=median(coverage),
        coverage.mean=mean(coverage),
        coverage.q75=quantile(coverage, 0.75),
        coverage.max=max(coverage),
        coverage.total=sum(coverage),
        bins.n.total=n(),
        bins.n.nz=sum(coverage > 0),
        bins.n.covered=sum(coverage > 1000)
    ) %>%
    ungroup() %>%
    mutate(
        bins.pct.nz=bins.n.nz / bins.n.total,
        bins.pct.covered=bins.n.covered / bins.n.total
    )
}

compute_coverage_summary <- function(
    filepath,
    ...){
    # Load bin-wise Interaction Frequency data
    coverage.df <- 
        filepath %>%
        load_genome_coverage()
    # compute summary stats of IF across the entire genome
    genome.df <- 
        coverage.df %>% 
        group_by(metric) %>% 
        compute_summary_stats() %>% 
        add_column(chr='Genome.Wide')
    # compute summary stats of IF per chr
    chrs.df <- 
        coverage.df %>% 
        group_by(chr, metric) %>% 
        compute_summary_stats()
    # Bind + melt summaries
    bind_rows(
        chrs.df, 
        genome.df
    ) %>% 
    dplyr::rename('count.type'=metric) %>% 
    pivot_longer(
        matches('^(coverage|bins)\\.'),
        names_to='metric',
        values_to='value'
    )
}

compute_all_coverage_summaries <- function(
    file.suffix='-coverage.tsv',
    resolutions=NULL,
    ...){
    # List input files (generated by cooltools coverage)
    parse_results_filelist(
        input_dir=COVERAGE_DIR,
        suffix=file.suffix,
        filename.column.name='MatrixID',
        param_delim='_'
    ) %>%
    mutate(MatrixID=str_remove(MatrixID, file.suffix)) %>% 
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
    # load each coverage file and count % of bins > 0 and > 1000 contacts
    mutate(
        coverage.summary=
            pmap(
                .l=.,
                .f=compute_coverage_summary,
                .progress=TRUE
            )
    ) %>%
    unnest(coverage.summary) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(count.type == 'cis') %>% 
    select(-c(filepath))
}

###################################################
# Plot Matrix QC stuff
###################################################
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
    ...){
    figure <- 
        plot.df %>% 
        ggplot(
            aes(
                x=range.start,
                y=value,
                color=orientation
            )
        ) + 
        geom_line() +
        scale_x_log10() +
        labs(
            x='Pair Distance',
            y='Raw Contacts'
        )
    figure %>% 
    post_process_plot(...)
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

