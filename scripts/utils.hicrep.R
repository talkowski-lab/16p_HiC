library(tidyverse)
library(magrittr)
library(glue)
library(tictoc)
###############
# Load resutls
load_all_hicrep_results <- function(){
    # Load all files generated from ./scripts/run.hicrep.sh
    parse_results_filelist(
        input_dir=HICREP_DIR,
        suffix='-hicrep.txt',
        filename.column.name='file.pair',
        param_delim='_'
    ) %>%
    separate_wider_delim(
        file.pair,
        delim='-',
        names=
            c(
                'A',
                'B',
                NA
            )
    ) %>% 
    # Get sample info for each pair
    separate_wider_delim(
        c(A, B),
        delim=fixed('.'),
        names_sep='.',
        names=
            c(
                "Edit", 
                "Genotype",
                "SampleNumber",
                "Celltype",
                NA,
                NA,
                "ReadFilter",
                NA
            )
    ) %>%
    mutate(
        A.isMerged=ifelse(grepl('Merged', A.SampleNumber), 'Merged', 'Individual'),
        B.isMerged=ifelse(grepl('Merged', B.SampleNumber), 'Merged', 'Individual')
    ) %>% 
    # Group HiCRep comparisons by mis/matching sample attributes
    pivot_longer(
        c(starts_with('A.'), starts_with('B.')),
        names_to='Sample.Attribute',
        values_to='Value'
    ) %>% 
    separate_wider_delim(
        Sample.Attribute,
        delim=fixed('.'),
        names=c('Sample.Index', 'Sample.Attribute')
    ) %>%
    pivot_wider(
        names_from=Sample.Index,
        values_from=Value
    ) %>% 
    rowwise() %>% 
    mutate(
        SamplePairType=
            case_when(
                A == B ~ A,
                A != B ~ 
                    c(A, B) %>% 
                    sort() %>%  
                    paste(collapse=" vs ")
            )
    ) %>% 
    select(-c(A, B)) %>%
    pivot_wider(
        names_from=Sample.Attribute,
        values_from=SamplePairType
    ) %>% 
    # Load hicrep scores for each comparison
    mutate(
        hicrep.results=
            purrr::pmap(
                .l=.,
                function(filepath, ...){
                    read_tsv(
                        filepath,
                        skip=2,
                        show_col_types=FALSE,
                        col_names=c('hicrep.score')
                    ) %>% 
                    add_column(chr=factor(CHROMOSOMES, levels=CHROMOSOMES))
                }
            )
    ) %>%
    unnest(hicrep.results) %>% 
    select(-c(filepath))
}
###############
# Plot resutls
plot_hicrep_boxplot <- function(
    hicrep_results,
    sample_group='Celltype',
    fill_var='is.downsampled', 
    facet_row=NULL,
    facet_col=NULL,
    scales='fixed',
    mark_df=NULL,
    mark_fill_var='chr',
    mark_alpha=0.6,
    mark_size=1,
    yintercept=0.95,
    ...){
    figure <- 
        hicrep_results %>%
        ggplot(
            aes(
                x=.data[[sample_group]],
                y=hicrep.score,
                fill=.data[[fill_var]]
            )
        ) +
        geom_hline(
            yintercept=yintercept,
            color='black',
            linetype='solid',
            linewidth=0.15
        ) +
        scale_y_continuous(labels=function(x) format(x, digits=2)) +
        theme(
            legend.position='top',
            axis.text.x=element_text(angle=35, hjust=1)
        ) +
        add_ggtheme()
    # Mark specific chromosomes in th boxplot if specified
    if (!(is.null(mark_df))) {
        figure <- 
            figure +
            geom_jitter(
                data=mark_df,
                aes(
                    x=.data[[sample_group]],
                    y=hicrep.score,
                    color=.data[[mark_fill_var]]
                ),
                alpha=mark_alpha,
                size=mark_size
            ) +
            geom_boxplot(outliers=FALSE)
    } else {
        figure <- 
            figure +
            geom_boxplot(outlier.size=0.5)
    }
    # Facet as specified
    if (!is.null(facet_col) & !is.null(facet_col)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                cols=vars(!!sym(facet_col)),
                scales=scales
            )
    } else if (!is.null(facet_row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                scales=scales
            )
    } else if (!is.null(facet_col)) {
        figure <- 
            figure +
            facet_grid2(
                cols=vars(!!sym(facet_col)),
                scales=scales
            )
    }
    figure
}
