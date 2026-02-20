###################################################
# Depdendencies
###################################################
library(here)
here::i_am('scripts/weiner.replication/make.replication.figures.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts',  'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR, 'scripts',  'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'weiner.replication', 'utils.weiner.replication.R'))
    library(magrittr)
    library(tidyverse)
})

###################################################
# Comparison parameters
###################################################
options(dplyr.summarize.inform=FALSE, scipen=999)
# parsed.args <- 
#     handle_CLI_args(
#         args=c('threads', 'force', 'resolutions'),
#         has.positional=FALSE
#     )
# Comparison hyper-params
resolutions <- c(100, 50, 10) * 1e3
hyper.params.df <- 
    expand_grid(
        # test.method=c('t.test', 'wilcox.test', 'anova', 'kruskal.test'),
        test.method=c('t.test', 'wilcox.test'),
        normalization=c('raw', 'balanced'),
        resolution=resolutions,
    ) %>% 
    cross_join(
        tribble(
             ~stat.var, ~facet.row, ~comparison.group, ~label.x.npc, ~label.y.npc, ~n_size, ~exclude_background,  ~ylim, ~width, ~height,  
               'group', 'Genotype',         'regions',      c(0.25),      c(0.48),       4,               FALSE,   TRUE,      8,       8,   
            'Genotype',    'group',        'Genotype',    c('left'),     c('top'),       3,                TRUE,  FALSE,      8,       5         
        )
    )
# list of regions to compare
region.comparisons <- 
    tribble(
        ~region.P2,    ~region.P1,     ~region.background,
        "16p11.2 CNV", "16p Telomere", "16p",
        "16p11.2 CNV", "16p Telomere", "16q",
        "16p11.2 CNV", "16p Telomere", "chr16",
        "16p11.2",     "16p Telomere", "16p",
        "16p11.2",     "16p Telomere", "16q",
        "16p11.2",     "16p Telomere", "chr16",
        # "16p",         "16p",          "16q",
        # "16q",         "16p",          "16p"
    )

###################################################
# Set up comparison inputs
###################################################
# list all comparisons across Samples, region-pairs and hyper-params
all.comparisons.df <- 
    list_all_region_comparisons(
        huper.params.df,
        region.comparisons
    ) %>% 
    mutate(
        resolution.label=scale_numbers(resolution, force_chr=TRUE),
        region.comparison=glue('{region.P2} - {region.P1} [{region.background}]'),
        results_file=
            file.path(
                ROBINSON_REPLICATION_DIR,
                'plots',
                glue('normalization_{normalization}'),
                glue('isMerged_{isMerged}'),
                glue('resolution_{resolution}'),
                glue('comparison_{comparison.group}'),
                glue('{region.P1}-{region.P2}-{region.background}-region.vs.distanct.matched.ctrl.pdf')
            ) %>%
            str_replace_all(' ', '.')
    ) %>% 
    # filter(test.method == 't.test') %>% 
    # filter(normalization == 'raw') %>% 
    filter(isMerged) %>% 
    filter(resolution %in% c(100000)) %>% 
    nest(
        plot.df=
            -c(
                stat.var, facet.row, comparison.group,
                results_file,
                region.comparison,
                resolution, resolution.label,
                normalization,
                test.method,
                isMerged
            )
    )

###################################################
# Compare regions differences
###################################################
all.comparisons.df %>%
    # filter(comparison.group == 'regions') %>% 
    pmap(
        .l=.,
        .f=plot_contacts_regions_boxplot,
        y.var='IF',
        facet.col='Celltype',
        # scales='free_y',
        # width=8,
        # height=8,
        # label.x.npc=c(0.25),
        # label.y.npc=c(0.48),
        # n_size=4,
        # exclude_background=FALSE,
        # ylim=TRUE,
        .progress=TRUE
    )

###################################################
# Compare genotype differences
###################################################
# all.comparisons.df %>%
#     filter(comparison.group == 'Genotype') %>% 
#     pmap(
#         .l=.,
#         .f=plot_contacts_regions_boxplot,
#         y.var='IF',
#         facet.col='Celltype',
#         scales='free_y',
#         # width=8,
#         # height=5,
#         # n_size=3,
#         # exclude_background=TRUE,
#         # ylim=FALSE,
#         .progress=TRUE
#     )
