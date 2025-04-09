library(glue)
library(tidyverse)
library(furrr)
library(multiHiCcompare)
library(viridis)
library(cowplot)
library(gtable)
# Dirs
BASE_DIR="/data/talkowski/Samples/16p_HiC"
source(glue('{BASE_DIR}/scripts/locations.R'))
source(glue('{BASE_DIR}/scripts/utils.data.R'))
PLOT_DIR=glue('{MULTIHICCOMPARE_DIR}/plots')
# Load sample info
sample_df <- 
    tribble(
        ~sample_ID, ~ celltype, ~gene, ~ genotype,
        '16pDELA3NSCHiC_S1.hg38.mapq_30.1000',   'NSC', '16p11', 'DEL',
        '16pDELB8NSCHiC_S2.hg38.mapq_30.1000',   'NSC', '16p11', 'DEL',
        '16pDUPC5NSCHiC_S3.hg38.mapq_30.1000',   'NSC', '16p11', 'DUP',
        '16pDUPD12NSCHiC_S4.hg38.mapq_30.1000',  'NSC', '16p11', 'DUP',
        '16pWTFACS1NSCHiC_S5.hg38.mapq_30.1000', 'NSC', '16p11', 'WT',
        '16pWTp46NSCHiC_S6.hg38.mapq_30.1000',   'NSC', '16p11', 'WT'
    ) %>%
    rowwise() %>% 
    mutate(
        sample_name=str_split(sample_ID, '\\.')[[1]][1],
        sample_group=paste(gene, genotype, sep='.')
    ) %>% 
    ungroup()
# All parameter combinations to run multHiCCompare for
param_combos <- 
    expand_grid(
        region=c(paste0('chr', 1:22), 'chrX'),
        resolution=c(500000, 100000, 50000, 10000, 5000),
        sample_groups=
            list(
              c('16p11.DEL', '16p11.WT'),
              c('16p11.DUP', '16p11.WT'),
              c('16p11.DUP', '16p11.DEL')
            ),
        zero.p=c(0, 0.2, 0.8, 1),
        A.min=c(0, 1, 5, 10),
        filter_bins=c(TRUE)
    ) %>% 
    rowwise() %>% 
    mutate(
        combo_name=
            glue("{paste(sample_groups, collapse='.vs.')}_{region}_{format(resolution, scientific=FALSE)}_{zero.p}_{A.min}_{filter_bins}"),
        output_file=glue("{MULTIHICCOMPARE_DIR}/results/{combo_name}.tsv"),
        md_plot_file=glue("{MULTIHICCOMPARE_DIR}/MD_plots/{combo_name}.pdf")
    ) %>% 
    ungroup()
param_combos %>% dim()
param_combos %>% select(-c(output_file, md_plot_file)) %>% head(2) %>% t()
# Run pipeline over param combos
# Parallel params
library(BiocParallel)
numCores <- 4
register(MulticoreParam(workers=numCores), default=TRUE)
data('hg38_cyto') # GRanges object with Centro/Telomere regions to filter
param_combos %>% 
filter(region == 'chrX') %>%
    filter(resolution == 10000) %>% 
    filter(grepl('16p11.*.vs.16p11.WT', combo_name)) %>% 
    filter(zero.p == 0.8) %>% 
    filter(A.min == 5) %>% 
    filter(filter_bins == TRUE) %>% 
    pmap(
         .l=.,
         .f=run_multiHiCCompare,
         sample_df=sample_df,
         .progress=TRUE
    )
# Load results all generated results
all.results <- 
    check_cached_results(
        results_file=glue('{MULTIHICCOMPARE_DIR}/16p.100kb.multiHiCCompare_results.tsv'),
        force_redo=FALSE,
        return_data=TRUE,
        results_fnc=load_all_multiHiCCompare_results,
        multiHiCCompare_results_dir=glue('{MULTIHICCOMPARE_DIR}/results'),
        resolution=100000,
        comparison='16p11.*.vs.16p11.WT',
        fdr.threshold=0.1
        # nom.threshold=0.05
    )
# Print results
all.results %>% head(2) %>% t()
# all.results %>% dplyr::count(zero.p, A.min, Is.Filtered, Resolution, Comparison)
# all.results %>% dplyr::count(Chr, chr, sig.lvl, Comparison)
# Clean data
plot_df <- 
    all.results %>% 
    mutate(logFC=-logFC) %>%
    mutate(log.p.adj=-log10(p.adj)) %>%
    left_join(load_chr_sizes(), by='Chr') %>% 
    rename(region1='region1.bp', region2='region2.bp') %>% 
    mutate(
        region1.bin=region1.bp / Resolution,
        region2.bin=region2.bp / Resolution,
        region1.pct=100 * (region1.bp / chr.total.bp),
        region2.pct=100 * (region2.bp / chr.total.bp)
    ) %>% 
    mutate(
        distance.bp=region2.bp - region1.bp,
        distance.kb=distance.bp / 1000,
        distance.bin=region2.bin - region1.bin,
        distance.pct=region2.pct - region1.pct
    ) %>%
    mutate(
        distance.discrete=
            case_when(
                distance.kb <= 5     ~  '    .<  5Kb',
                distance.kb <= 20    ~ ' 5Kb<.<20Kb',
                distance.kb <= 50    ~ '20Kb<.<50Kb',
                distance.kb <= 1000  ~ '50Kb<.< 1Mb',
                distance.kb <= 5000  ~ ' 1Mb<.<50Kb',
                distance.kb <= 30000 ~ ' 5Mb<.<30Kb',
                distance.kb >  30000 ~ '30Mb<.     ',
            ) %>%
            factor(levels=c('    .<  5Kb', ' 5Kb<.<20Kb', '20Kb<.<50Kb', '50Kb<.< 1Mb', ' 1Mb<.<50Kb', ' 5Mb<.<30Kb', '30Mb<.     '))
    ) %>% 
    mutate(
        Chr=factor(Chr, levels=c(paste0('chr', 1:22), 'chrX')),
        Resolution=fct_reorder(glue('{Resolution / 1000}Kb'), Resolution)
    )
plot_df %>% head(2) %>% t()
# Plots
# Regular Manhattan plot
manhattan_regular <- 
    plot_df %>% 
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    mutate(tmp=as.integer(Chr) %% 2 == 0) %>% 
    ggplot(
        aes(
           x=Chr,
           y=-log10(p.adj),
           color=Comparison,
        ),
        size=0.8,
        alpha=0.3
    ) +
    geom_jitter(size=1) +
    geom_vline(
        aes(
            xintercept=as.integer(Chr) - 0.5,
        ),
        linewidth=0.1,
        linetype='dashed'
    ) +
    # geom_vline(aes(xintercept=stage(x, after_scale = xintercept - 0.5)), linetype='dashed') +
    labs(title='Genome-Wide differential contacts from 16p') +
    scale_y_continuous(
        expand=c(0.01, 0.01, 0.01, 0.01),
        limits=c(0, NA)
    ) +
    theme(axis.title.x=element_blank(), legend.position='top') +
    add_ggtheme()
ggsave(
    glue('{PLOT_DIR}/genome.regular.manhattan_plot.pdf'),
    manhattan_regular, width=10, height=6, units='in'
)
## Volcano Plots genome wide
plot_df %>%
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    multiHiCCompare_genome_volcano_plot(
        output_file=glue('{PLOT_DIR}/volcano/genome.volcano_plot.pdf'),
        color='distance.kb',
        scales='fixed',
        pal.direction=-1,
        pal='A',
        ncol=1,
        size=1,
        width=8,
        height=7
    )
## Volcano Plots all chromosomes
plot_df %>%
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    multiHiCCompare_allChromosome_volcano_plot(
        output_file=glue('{PLOT_DIR}/volcano/all.chromosomes.volcano_plot.pdf'),
        color='distance.discrete',
        pal.direction=-1,
        pal='A',
        ncol=3,
        size=1,
        width=14,
        height=20
    )
## Volcano Plots per chromosome
plot_df %>%
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    group_split(chr) %>%
    purrr::map(~  
        multiHiCCompare_chromosome_volcano_plot(
            .,
            chr=.$chr[[1]],
            output_file=glue('{PLOT_DIR}/volcano/chr{.$chr[[1]]}.volcano_plot.pdf'),
            width=7,
            height=7,
        )
    )
## Manhattan-esque plot all chromosomes
plot_df %>%
    pivot_longer(
        c(logCPM, logFC, log.p.adj), 
        names_to='statistic',
        values_to='value',
    ) %>%
    select(-c(logCPM)) %>% 
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    group_split(chr, Resolution) %>%
    purrr::map(~  
        multiHiCCompare_manhattan_plot(
            .,
            chr=.$chr[[1]],
            Resolution=.$Resolution[[1]],
            y_axis='region1.band',
            output_file=glue('{PLOT_DIR}/manhattan/all.chrs.{.$Resolution[[1]]}.manhattan.pdf'),
            width=7,
            height=15,
        )
    )
## Manhattan-esque plot per chromosome
plot_df %>%
    pivot_longer(
        c(logCPM, logFC, log.p.adj), 
        names_to='statistic',
        values_to='value',
    ) %>%
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    group_split(chr, Resolution) %>%
    purrr::map(~  
        multiHiCCompare_manhattan_plot(
            .,
            chr=.$chr[[1]],
            Resolution=.$Resolution[[1]],
            y_axis='region1.bin',
            output_file=glue('{PLOT_DIR}/manhattan/chr{.$chr[[1]]}.{.$Resolution[[1]]}.manhattan.pdf'),
            width=7,
            height=15,
        )
    )
## genome distance vs FC scatter plot
plot_df %>%
    select(-c(distance.discrete)) %>% 
    pivot_longer(
        starts_with('distance.'),
        names_to='distance.unit',
        values_to='distance.value',
    ) %>%
    mutate(Comparison=str_remove(Comparison, '.vs.16p11.WT')) %>% 
    multiHiCCompare_genome_scatter_plot(
        output_file=glue('{PLOT_DIR}/scatter/genome.dist_scatter.pdf'),
        ncol=1,
        size=1,
        width=8,
        height=7
    )
