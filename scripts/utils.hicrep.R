###################################################
# Dependencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(tictoc)

###################################################
# Load resutls
###################################################
load_all_hicrep_results <- function(sample_metadata=NULL){
    # Load all files generated from ./scripts/run.hicrep.sh
    # sample_metadata=sample.metadata %>% select( SampleID, Batch)
    parse_results_filelist(
        input_dir=HICREP_DIR,
        suffix='-hicrep.txt',
        filename.column.name='file.pair',
        param_delim='_'
    ) %>%
    # Separate IDs of 2 matrices being compared for each results file
    separate_wider_delim(
        file.pair,
        delim='-',
        names=
            c(
                'MatrixID.P1',
                'MatrixID.P2',
                NA
            )
    ) %>% 
    # Extract sample metadata from IDs
    get_info_from_MatrixIDs(
        matrix_ID_col='MatrixID.P1',
        keep_id=FALSE,
        sample_ID_col='SampleInfo.P1.SampleID',
        col_prefix='SampleInfo.P1.',
        nest_col=NA
    ) %>% 
    # Add extra sample metadata as paired columns 
    # NOTE: using right_join to filter samples only present in the metadata table
    inner_join(
        .,
        sample_metadata %>% 
        rename_with(
            .fn=~ str_replace(.x, '^', 'SampleInfo.P1.'),
            .cols=-c(SampleID)
        ),
        by=join_by(SampleInfo.P1.SampleID == SampleID)
    ) %>% 
    # Repeat for the other sample in each pair
    get_info_from_MatrixIDs(
        matrix_ID_col='MatrixID.P2',
        keep_id=FALSE,
        sample_ID_col='SampleInfo.P2.SampleID',
        col_prefix='SampleInfo.P2.',
        nest_col=NA
    ) %>% 
    # NOTE: using right_join to filter samples only present in the metadata table
    inner_join(
        .,
        sample_metadata %>% 
        rename_with(
            .fn=~ str_replace(.x, '^', 'SampleInfo.P2.'),
            .cols=-c(SampleID)
        ),
        by=join_by(SampleInfo.P2.SampleID == SampleID)
    ) %>% 
    # Now format sample metadata per pair for easy grouping+plotting
    nest(
        SampleInfo.P1=starts_with('SampleInfo.P1.'),
        SampleInfo.P2=starts_with('SampleInfo.P2.')
    ) %>%
    rowwise() %>% 
    mutate(
        SamplePairInfo=
            merge_sample_info(
                SampleInfo.P1,
                SampleInfo.P2,
                prefix.P1='SampleInfo.P1.',
                prefix.P2='SampleInfo.P2.'
            ) %>%
            list()
    ) %>% 
    ungroup() %>% 
    select(-c(starts_with('SampleInfo'))) %>% 
    unnest(SamplePairInfo) %>% 
    # Annotate which resolution is the minimum viable for each SampleID
    # get_min_resolution_per_matrix(filter_res=FALSE) %>% 
    # # look if resolution is "ideal" according to original HiCRep authors
    # left_join(
    #     RESOLUTION_IDEAL_H,
    #     relationship='many-to-many',
    #     by=join_by(resolution)
    # ) %>% 
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
                    # every line is a score per chromosome in order
                    add_column(chr=factor(CHROMOSOMES, levels=CHROMOSOMES))
                },
                .progress=TRUE
            )
    ) %>%
    unnest(hicrep.results)
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

###################################################
# Plot resutls
###################################################
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
