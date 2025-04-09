library(glue)
library(tidyverse)
library(furrr)
library(ggh4x)
library(hictkR)
library(ggnewscale)
# library(parallel)
###############
# Dirs
# BASE_DIR="/data/talkowski/Samples/16p_HiC"
BASE_DIR="/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p"
source(glue('{BASE_DIR}/scripts/utils.data.R'))
PAIRS_DIR=glue('{BASE_DIR}/results/mapped_parsed_sorted_chunks')
RESULTS_DIR=glue('{MATRIX_QC_DIR}/bin.summaries')
PLOT_DIR=glue('{MATRIX_QC_DIR}/plots')
###############
# Params
RESOLUTIONS=c(10000, 50000, 100000, 500000, 1000000)
RESOLUTION_IDEAL_h <- 
    c(
        '10000'=20,
        '25000'=10,
        '40000'=5,
        '100000'=3,
        '500000'=2,
        '500000'=1,
        '1000000'=1,
        '1000000'=0
    )
NORMALIZATIONS=c('cis.ice', 'NONE')
###############
# Load pairtools stats
sample.metadata <- load_sample_metadata()
pair.stats.list <- load_pairtools_stats(pairs_dir=PAIRS_DIR)
###############
# Plot distance by by frequency for QC reasons
plot.df <- 
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
                orientation =='-+'   & range.end < 1001  ~ 'self-circle', 
                orientation =='+-'   & range.end < 1001  ~ 'dangling-end', 
                range.start >= 20000                     ~ 'Long.Range  >  20Kb',
                range.start >=  1000 & range.start < 20001 ~ 'Short.Range <= 20Kb',
                range.end   <=  1000                     ~ 'Too.Short   <=  1Kb',
                interaction == 'cis' & range.end == Inf ~ 'Cis.Total',
                interaction == 'trans'                   ~ 'Trans',
            )
    ) %>%
    mutate(
        Category=
            factor(
                Category,
                levels=
                    c(
                        'self-circle', 
                        'dangling-end', 
                        'Too.Short   <=  1Kb',
                        'Short.Range <= 20Kb',
                        'Long.Range  >  20Kb',
                        'Cis.Total',
                        'Trans'
                    )
            )
    ) %>% 
    group_by(SampleID, Category) %>% 
    summarize(
        value=sum(value),
        # total per sample
        total.unique.contacts=unique(total.unique.contacts)
    ) %>% 
    mutate(frac=100 * (value / total.unique.contacts))
plot.df  %>%
    plot_distance_freqs(y_val='frac') %>%
    ggsave(
        glue('{PLOT_DIR}/distance.contact.frequency.barplot.pdf'),
        .,
        height=6,
        width=10,
        units='in'
    )
# Plot chrom level heatmap
pair.stats.list$chr.stats %>% 
    mutate(
        pair.types=
            case_when( 
                grepl('(_|EBV)', chr1) | grepl('(_|EBV)', chr2) ~ 'artifact',
                chr1 == chr2 ~ 'cis',
                chr1 != chr2 ~ 'trans'
            )
    ) %>%
    group_by(SampleID, pair.types) %>% 
    dplyr::summarize(sum=sum(value))
plot.df  <- 
    pair.stats.list$chr.stats %>% 
    filter( 
        if_all(
            starts_with('chr'),
            ~ !grepl('(_|EBV)', .x)
        )
    ) %>% 
    mutate(across(starts_with('chr'), ~ factor(.x, levels=CHROMOSOMES))) %>% 
    rowwise() %>%
    mutate(
        `% Total Contacts`=100 * (value / total.unique.contacts),
        `log10(Contacts)`=log10(value)
    ) #%>%
    # pivot_longer(c(contacts, pct), names_to='metric', values_to='value')
plot_heatmap <- function(plot.df, fill_var){
    plot.df %>%
    ggplot(aes(x=chr1, y=chr2,)) +
    geom_tile(aes(fill=.data[[fill_var]])) +
    scale_fill_gradient(low='grey90', high='red') +
    facet_wrap(~ SampleID) +
    theme(legend.position='top', axis.text.x=element_text(angle=45, hjust=1)) +
    add_ggtheme()
}
plot.df %>% 
    plot_heatmap(fill_var="log10(Contacts)") %>% 
    ggsave(
        glue('{PLOT_DIR}/chrom.Inter.Intra.log10.heatmap.pdf'),
        .,
        height=6,
        width=10,
        units='in'
    )
plot.df %>% 
    plot_heatmap(fill_var="% Total Contacts") %>% 
    ggsave(
        glue('{PLOT_DIR}/chrom.Inter.Intra.pct.heatmap.pdf'),
        .,
        height=6,
        width=10,
        units='in'
    )
###############
# Load HiC Matrices, get per-bin summary stats
bin.summaries <- 
    load_sample_metadata() %>% 
    add_column(dummy=1) %>% 
    full_join(
        expand_grid(
            normalization=NORMALIZATIONS,
            resolution=RESOLUTIONS,
            dummy=1
        ),
        by='dummy',
        relationship='many-to-many',
    ) %>%  
    mutate(results_file=glue('{RESULTS_DIR}/{format(resolution, trim=TRUE, scientific=FALSE)}/{normalization}/{Sample.Matrix}-bin_summaries.tsv')) %>% 
    select(-c(dummy)) %>% 
    arrange(desc(resolution), desc(normalization), Sample.ID) %>%
    # filter(resolution %in% c('gw.ice')) %>% 
    # filter(resolution < 50000) %>% 
    mutate(
        bin.summaries=
            pmap(
                .l=.,
                .f=check_cached_results,
                force_redo=FALSE,
                return_data=TRUE,
                show_col_types=FALSE,
                results_fnc=summarize_contacts_per_bin,
                cluster=NULL,
                .progress=TRUE
            )
    ) %>%
    unnest(bin.summaries) %>% 
    pivot_wider(names_from=stat, values_from=value) %>% 
    select(-c(results_file, mcool_file))
# bin.summaries %>% dplyr::count(SampleID, ReadFilter, normalization, resolution, Chr)
bin.summaries %>% head(2) %>% t()
plot.df <- 
    bin.summaries %>%
    select(-c(ReadFilter, SampleMatrix)) %>% 
    # filter(resolution == 10000) %>%
    # filter(normalization == 'NONE') %>% 
    # filter(SampleID == "16pDELA3NSCHiC_S1") %>% 
    mutate(Chr=factor(Chr, levels=CHROMOSOMES)) %>% 
    group_by(normalization, resolution, SampleID) %>% 
    arrange(Chr, bin.start) %>%
    mutate(bin.id=row_number()) %>%
    ungroup() %>%
    mutate(resolution=glue('{format(resolution / 1000, scientific=FALSE, trim=TRUE)}Kb')) %>% 
    mutate(resolution=factor(resolution, levels=RESOLUTIONS))
    mutate(normalization=case_when(normalization == 'NONE' ~ 'Raw', normalization == 'cis.ice' ~ 'ICE'))
plot.df %>% arrange(SampleID, Chr, bin.start) %>% head(2) %>% t()
plot.df %>% 
    filter(!is.merged) %>% 
    plot_bin_totals_stats() %>%
    ggsave(glue('{PLOT_DIR}/genome.bin_totals.boxplot.pdf'), ., width=10, height=10, units='in')
plot.df %>%
    filter(!is.merged) %>% 
    filter(normalization != 'cis.vc') %>% 
    filter(Chr == 'chr16') %>% 
    # filter(bin.start < 38e6 & bin.start > 20e6) %>% 
    filter(bin.start < 38e6 & bin.start > 0) %>% 
    plot_totals_across_chr16() %>% 
    ggsave(glue('{PLOT_DIR}/chr16.lineplot.pdf'), ., width=10, height=7, units='in')
plot.df %>%
    filter(!is.merged) %>% 
    filter(normalization != 'cis.vc') %>% 
    plot_bin_stats_across_region() %>% 
    ggsave(glue('{PLOT_DIR}/genome.lineplot.pdf'), ., width=10, height=10, units='in')
###############
# Plot coverages genome-wide @ diff resolutions (bins - coverage)
## Plot with/without mapq_30 filtering
## Plot with/without filtering
## Plot before/after balancing
###############
# Load hicrep data
hicrep.results <- 
    check_cached_results(
        results_file=glue('{MATRIX_QC_DIR}/hicrep.results.tsv'),
        force_redo=TRUE,
        return_data=TRUE,
        results_fnc=load_hicrep_results,
        resolutions=NULL
        # resolutions=c(10000, 100000)
    )
hicrep_results %>% head(2) %>% t()
hicrep_results %>% 
    # filter(Chr == 'chr1') %>% 
    # filter(ReadFilter == 'mapq_30') %>% 
    # filter(Resolution == '100000') %>% 
    # filter(Genotype.Pair == 'DEL.vs.DEL') %>% 
    # dplyr::count(is.Downsampled, h, Max.Window.Size) %>%
    # dplyr::count(A.SampleID, B.SampleID, ReadFilter)
    # distinct(ReadFilter, Resolution, is.Downsampled, h, Max.Window.Size) %>% 
    distinct(ReadFilter, Resolution, h, Max.Window.Size) %>% 
    # dplyr::count(ReadFilter, Resolution, is.Downsampled, Genotype.Pair) %>%
    # dplyr::count(is.Downsampled, Genotype.Pair) %>%
    print(n=Inf)
is_h_ideal <- function(Resolution, h){
    if (Resolution %in% names(RESOLUTION_IDEAL_h)) {
        if (RESOLUTION_IDEAL_h[[Resolution]] == h) {
            return('Ideal')
        } else {
            return('Non.Ideal')
        }
    } else {
        return('No Ideal H')
    }
}
plot_df <- 
    hicrep_results %>%
    filter(ReadFilter != 'mixed') %>% 
    filter(!(grepl('_S[0-9]S[0-9]', A.SampleID) | grepl('_S[0-9]S[0-9]', B.SampleID))) %>%
    rowwise() %>% 
    mutate(h.ideal=is_h_ideal(Resolution, h)) %>%
    ungroup()
plot_hicrep_results <- function(hicrep_results, size=0.7){
    hicrep_results %>%
    ggplot(
        aes(
            x=Genotype.Pair,
            y=hicrep.score,
            fill=is.Downsampled
        )
    ) +
    geom_boxplot(outlier.size=0.5) +
    geom_jitter(
        data=. %>% filter(h.ideal == 'Ideal'),
        aes(
            x=Genotype.Pair,
            y=hicrep.score,
            # color=h.ideal,
            fill=is.Downsampled
        ),
        color='green',
        size=size
    ) +
    facet_grid2(rows=vars(Resolution), cols=vars(ReadFilter), scales='fixed') +
    theme(
        legend.position='top', 
        axis.text.x=element_text(angle=45, hjust=1)
    ) +
    add_ggtheme()
}
plot_df %>% 
    plot_hicrep_results() %>%
    ggsave(glue('{PLOT_DIR}/hicrep.pdf'), ., width=10, height=10, units='in')
plot_df %>% 
    filter(Chr == 'chr16') %>% 
    plot_hicrep_results(size=2) %>%
    ggsave(glue('{PLOT_DIR}/hicrep.chr16.pdf'), ., width=10, height=10, units='in')
plot_df %>% 
    filter(Chr == 'chr9') %>% 
    plot_hicrep_results(size=2) %>%
    ggsave(glue('{PLOT_DIR}/hicrep.chr9.pdf'), ., width=10, height=10, units='in')
