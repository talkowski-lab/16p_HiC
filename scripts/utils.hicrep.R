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
            ),
        cols_remove=FALSE
    ) %>%
    rename(
        'A.SampleID'=A.A,
        'B.SampleID'=B.B
    ) %>% 
    mutate(
        across(
            ends_with('SampleID'),
             ~ .x %>%
               str_extract('(([^\\.]+\\.){3}[^\\.]+)\\.(.*$)', group=1) %>%
               str_replace_all(fixed('.'), '~') # temp change delim so ID stays intact
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
    # Now group pairs by difference in specific sample metadata (genotype, celltype)
    rowwise() %>% 
    mutate(
        SamplePairType=
            case_when(
                A == B ~ glue('{A} vs {A}'),
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
    # Fix Sample.ID delim
    mutate(SampleID=str_replace_all(SampleID, fixed('~'), fixed('.'))) %>% 
    rename('Sample.ID'=SampleID) %>% 
    # Load hicrep scores for each comparison
    mutate(
        hicrep.results=
            purrr::pmap(
                .l=.,
                function(filepath, ...){
                    read_tsv(
                        filepath,
                        skip=2,
                        progress=FALSE,
                        show_col_types=FALSE,
                        col_names=c('hicrep.score')
                    ) %>% 
                    add_column(chr=factor(CHROMOSOMES, levels=CHROMOSOMES))
                },
                .progress=TRUE
            )
    ) %>%
    unnest(hicrep.results) %>% 
    select(-c(filepath))
}

make_heatmap_plotdf <- function(
    hicrep.results,
    ...){
    hicrep.results %>% 
    group_by(resolution, Sample.ID) %>% 
    summarize(`Mean HiCRep Score`=mean(hicrep.score)) %>%
    ungroup() %>% 
    # Get sample metadata for both samples per pair
    separate_wider_delim(
        Sample.ID,
        delim=' vs ',
        names=c('A.Sample.ID', 'B.Sample.ID')
    ) %>%
    separate_wider_delim(
        c(A.Sample.ID, B.Sample.ID),
        delim=fixed('.'),
        names_sep='.',
        names=c('Edit', 'Genotype', 'SampleNumber', 'Celltype'),
        cols_remove=FALSE
    ) %>%
    rename(
        'A.Sample.ID'=A.Sample.ID.A.Sample.ID,
        'B.Sample.ID'=B.Sample.ID.B.Sample.ID
    ) %>% 
    rename_with(
        .fn=~ str_replace(.x, '(A|B)\\.Sample\\.ID\\.(.*)$', '\\1.\\2'),
        .cols=matches('^(A|B).Sample.ID.')
    )
}
###############
# Plot resutls
plot_hicrep_boxplot <- function(
    plot.df,
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
    # Make base boxplot
    if (is.na(fill_var)) {
        figure <- 
            ggplot(
                plot.df,
                aes(
                    x=.data[[sample_group]],
                    y=hicrep.score
                )
            )
    } else {
        figure <- 
            ggplot(
                plot.df,
                aes(
                    x=.data[[sample_group]],
                    y=hicrep.score,
                    fill=.data[[fill_var]]
                )
            )
    }
    # Mark specific chromosomes in th boxplot if specified
    if (!(is.null(mark_df))) {
        figure <- 
            figure +
            geom_boxplot(outliers=FALSE) +
            geom_jitter(
                data=mark_df,
                aes(
                    x=.data[[sample_group]],
                    y=hicrep.score,
                    color=.data[[mark_fill_var]]
                ),
                alpha=mark_alpha,
                size=mark_size
            )
    } else {
        figure <- 
            figure +
            geom_boxplot(outlier.size=0.5)
    }
    # Add reference line for comparison
    figure <- 
        figure + 
        geom_hline(
            yintercept=yintercept,
            color='black',
            linetype='solid',
            linewidth=0.15
        )
        # scale_y_continuous(labels=function(x) format(x, digits=2)) +
    # Add faceting + theming
    figure <- 
        figure %>% 
        post_process_plot(
            facet_col=facet_row,
            facet_row=facet_col,
            scales=scales,
            scale_mode='',
            log_base=10,
            axis_label_accuracy=2,
            legend.position='top',
            axis.text.x=element_text(angle=35, hjust=1),
            ...
            # n_breaks=NULL,
            # limits=NULL,
            # expand=c(0.00, 0.00, 0.00, 0.00),
            # legend.position='right',
            # legend.ncols=1,
            # axis.text.x=element_text(angle=35, hjust=1),
        )
}

plot_hicrep_heatmap <- function(
    plot.df,
    x_var='A.chr',
    y_var='B.chr',
    fill_var='Mean.HiCRep.Score',
    x_text_angle=25,
    ...){
    figure <- 
        plot.df %>%
        ggplot(
            aes(
                x=.data[[x_var]],
                y=.data[[y_var]]
            )
        ) +
        geom_tile(
            aes(
                fill=.data[[fill_var]]
            )
        ) +
        scale_fill_gradient(
            low='grey90',
            high='blue'
        ) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        facet_nested(
            rows=vars(B.Edit, B.Genotype),
            cols=vars(A.Edit, A.Genotype),
            scales='free',
            space='free',
            nest_line=
                element_line(
                    color='black',
                    linetype='solid',
                    linewidth=0.40
                ),
        ) +
        make_ggtheme(
            ...,
            legend.position='right',
            panel.spacing=unit(0, "lines"),
            panel.border=
                element_rect(
                    fill=NA,
                    color='black',
                    linetype='solid',
                    linewidth=0.40
                ),
            strip.background=
                element_rect(
                    color='black',
                    linetype='solid',
                    linewidth=0.40
                ),
            strip.text=element_text(face='bold', size=10),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(angle=x_text_angle, hjust=1)
        ) 
    figure
}
