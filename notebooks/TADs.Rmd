library(glue)
library(furrr)
library(tidyverse)
library(HiContacts)
library(viridis)
library(ggplot2)
library(ggh4x)
# library(ggtext)
library(cowplot)

# BASE_DIR="/data/talkowski/broadIncoming/22LCC2LT4/fastq"
# BASE_DIR="/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p"
BASE_DIR="/data/talkowski/Samples/16p_HiC"
PLOT_DIR=glue('{BASE_DIR}/results/TADs/plots')
source(glue('{BASE_DIR}/scripts/locations.R'))
source(glue('{BASE_DIR}/scripts/utils.data.R'))
# Load Directionality Index scores
DI.scores <- 
    load_all_TAD_annotations(
        tad_annotations_dir=glue('{BASE_DIR}/results/TADs/hiTAD'),
        tad.method='hiTAD-DIs',
        file_suffix='-DIs.txt',
        param_names=
            c(
                'Normalization',
                'Resolution',
                'SampleName'
            ),
    ) %>%
    separate_wider_delim(
        SampleName,
        delim=fixed('.'),
        names=
            c(
                'SampleName',
                NA,
                'ReadFilter',
                NA
            )
    ) %>% 
    mutate(
        Resolution=as.integer(Resolution),
        chr=factor(chr, levels=c(paste0('chr', 1:22), 'chrX')),
        Genotype=
            case_when(
                str_detect(SampleName, '^16pWT') ~ 'WT',
                str_detect(SampleName, '^16pDEL') ~ 'DEL',
                str_detect(SampleName, '^16pDUP') ~ 'DUP',
                TRUE ~ 'Unknown'
            ) %>%
            factor(levels=c('WT', 'DEL', 'DUP')),
        Replicate=
            SampleName %>%
            str_extract('^16p.*_S([1-6])\\.*', group=1) %>% 
            as.integer() %>% 
            { 2 - . %% 2 } %>% 
            as.factor()
    )
DI.scores %>% head(2) %>% t()
# Load all TAD annotations
tad.annotations <- 
    load_all_TAD_annotations(
        tad_annotations_dir=glue('{BASE_DIR}/results/TADs/hiTAD'),
        tad.method='hiTAD',
        file_suffix='-tads.txt',
        param_names=
            c(
                'Normalization',
                'Resolution',
                'SampleName'
            ),
    ) %>%
    separate_wider_delim(
        SampleName,
        delim=fixed('.'),
        names=
            c(
                'SampleName',
                NA,
                'ReadFilter',
                NA
            )
    ) %>% 
    # Get chromosome size in bins per resolution for convienience
    left_join(
        DI.scores %>%
        select(-c(DI, Method)) %>%
        group_by(across(-starts_with('bin.'))) %>%
        dplyr::count(name='chrom.total.bins')
    ) %>%
    mutate(
        Resolution=as.integer(Resolution),
        chr=factor(chr, levels=c(paste0('chr', 1:22), 'chrX')),
        Genotype=
            case_when(
                str_detect(SampleName, '16pWT') ~ 'WT',
                str_detect(SampleName, '16pDEL') ~ 'DEL',
                str_detect(SampleName, '16pDUP') ~ 'DUP',
                TRUE ~ 'Unknown'
            ) %>%
            factor(levels=c('WT', 'DEL', 'DUP')),
        Replicate=
            SampleName %>%
            str_extract('^16p.*_S([1-6])\\.*', group=1) %>% 
            as.integer() %>% 
            { 2 - . %% 2 } %>% 
            as.factor()
    ) 
tad.annotations %>% head(2) %>% t()
tad.annotations %>% dplyr::count(Method, Normalization, ReadFilter, Resolution, SampleName, Genotype) %>% print(n=Inf)
# Plot # of TADs across samples
tad.annotations %>%
    dplyr::count(Resolution, Genotype, SampleName, chr) %>%
    group_by(Resolution) %>%
    group_map(
        ~ ggplot(
            .x, 
            aes(
                x=chr,
                y=SampleName,
                fill=n
            )
        ) +
        geom_tile() +
        scale_fill_viridis(discrete=FALSE) +
        labs(title=sprintf('16p sample TADs @ %sKb', .x$Resolution / 1000)) +
        theme(
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_markdown()
        ) +
        add_ggtheme(),
        .keep=TRUE
    ) %>% 
    cowplot::plot_grid(
        plotlist=.,
        align='v',
        ncol=2
    ) %>% 
    ggsave(
        glue('{PLOT_DIR}/16p.hiTAD.TAD_Counts.pdf'),
        plot=.,
        height=6,
        width=10,
        units='in'
    )
# Plot TAD length Distribution
tad.annotations %>%
    filter(chr %in% paste0('chr', c(1,3,5,10,16,20,'X'))) %>% 
    mutate(TAD.length=(TAD.end - TAD.start) / Resolution) %>% 
    mutate(Resolution=glue('{Resolution / 1000}Kb')) %>% 
    mutate(Resolution=factor(Resolution, levels=c('100Kb', '50Kb', '10Kb', '5Kb'))) %>% 
    group_by(Method) %>%
    group_map(
        ~ ggplot(
            .x, 
            aes(
                y=TAD.length,
                x=Genotype,
                color=Replicate,
            )
        ) +
        geom_boxplot() +
        labs(x='TAD Length (bins)') +
        ggh4x::facet_grid2(cols=vars(chr), rows=vars(Resolution), scales='free_y') + #, independent='y') +
        theme(
            axis.title.x=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_markdown()
        ) +
        add_ggtheme(),
        .keep=TRUE
    ) %>% 
    cowplot::plot_grid(
        plotlist=.,
        align='v',
        ncol=1
    ) %>% 
    ggsave(
        glue('{PLOT_DIR}/16p.hiTAD.TAD_Lengths.pdf'),
        plot=.,
        height=10,
        width=20,
        units='in'
    )
# Plot DI Scores
tmp <- 
    DI.scores %>%
    add_count(Resolution, SampleName, chr, name='n.bins') %>% 
    left_join(
        tad.annotations %>% select(c(Resolution, SampleName, chr, TAD.start, TAD.end)),
        by=
            join_by(
                Resolution,
                chr,
                SampleName,
                between(bin.start, TAD.start, TAD.end, bounds='[)')
            ) 
    ) %>%
    mutate(
        TAD.boundaries=
            case_when(
                bin.start == TAD.start ~ 'TAD.Start',
                bin.end == TAD.end ~ 'TAD.End',
                bin.start >= TAD.start & bin.end <= TAD.end ~ 'TAD.inner',
                TRUE ~ 'Non-TAD'
        )
    )
tmp %>% 
    # filter(chr %in% paste0('chr', c(1,3,5,10,16,20,'X'))) %>% 
    filter(chr == 'chr16') %>% 
    mutate(Resolution=glue('{Resolution / 1000}Kb')) %>%
    mutate(Resolution=factor(Resolution, levels=c('100Kb', '50Kb', '10Kb', '5Kb'))) %>% 
    group_by(Method) %>%
    group_map(
        ~ ggplot(
            .x, 
            aes(
                x=bin.start,
                y=DI,
                group=chr,
                color=Genotype
                # color=Genotype, shape=Replicate
            )
        ) +
        # geom_rect(data=.x %>% filter(SampleName == ''), aes(xmin=TAD.start, xmax=TAD.end, ymin=-Inf, ymax=0, fill=SampleName), alpha=0.01, size=0) +
        geom_line() +
        labs(x='bins', y='Directionality Index') +
        ggh4x::facet_grid2(cols=vars(chr), rows=vars(Resolution), scales='free_y') +
        scale_x_continuous(
            n.breaks=10,
            expand=c(0.01, 0.01, 0.01, 0.01),
            limits=c(0, NA)
        ) +
        theme(
            legend.position='top',
            axis.title.x=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_markdown()
        ) +
        add_ggtheme(),
        .keep=TRUE
    ) %>% 
    cowplot::plot_grid(
        plotlist=.,
        align='v',
        ncol=1
    ) %>% 
    ggsave(
        glue('{PLOT_DIR}/16p.hiTAD.DI_values.pdf'),
        plot=.,
        height=10,
        width=20,
        units='in'
    )






cads.test <- 
    tad.annotations %>%
    select(-c(Normalization, ReadFilter, Method)) %>%
    filter(Resolution == '100000') %>%
    filter(chr == 'chr16')
tmp <- 
    DI.scores %>%
    add_count(Resolution, SampleName, chr, name='n.bins') %>% 
    filter(Resolution == '100000') %>%
    filter(chr == 'chr16') %>% 
    inner_join(
        tads.test,
        by=join_by(Resolution, chr, SampleName, between(bin.start, TAD.start, TAD.end, bounds='[)')) 
    ) %>%
    arrange(chr, bin.start)
tmp %>% group_by(SampleName) %>% summarize(n=n(), total=unique(n.bins)) 
tads.test %>% mutate(len=(TAD.end - TAD.start) / Resolution) %>% group_by(SampleName) %>% summarize(tad.covered=sum(len))
tmp %>% group_by(SampleName) %>% summarize(n=n(), total=unique(n.bins))
    # group_by(SampleName, chr, Resolution) %>% 
    # summarize(TAD.starts=list(TAD.start), TAD.ends=list(TAD.end)) %>%
    # mutate(TAD.boundaries=list(sort(c(unlist(TAD.starts), unlist(TAD.ends))))) %>%
    # mutate(nonTAD.boundaries=

