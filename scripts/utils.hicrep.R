library(tictoc)
library(glue)
# Match ideal smoothing param to specified resolution based on this
# https://github.com/TaoYang-dev/hicrep?tab=readme-ov-file#hicrep-parameters
HICREP_PARAM_NAMES <- 
    c(
        'resolution',
        'h',
        'is.downsampled',
        'window.size'
    )

load_all_hicrep_results <- function(
    input_dir=HICREP_DIR,
    param_names=HICREP_PARAM_NAMES,
    suffix='hicrep.txt',
    delim='-'){
    # Load all files generated from ./scripts/run.hicrep.shS
    input_dir %>% 
    list.files(
        recursive=TRUE,
        pattern=glue('*{suffix}'),
        full.names=FALSE
    ) %>%
    tibble(fileinfo=.) %>%
    mutate(filepath=file.path(input_dir, fileinfo)) %>%
    # Get hicrep params
    separate_wider_delim(
        fileinfo,
        delim='/',
        names=c(param_names, 'file.pair')
    ) %>% 
    rowwise() %>% 
    mutate(
        across(
            param_names,
            ~ str_remove(
                .x,
                pattern=glue('({paste(param_names, collapse="|")})_')
            )
        )
    ) %>% 
    separate_wider_delim(
        file.pair,
        delim=delim,
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
plot_hicrep_boxplot <- function(hicrep_results, size=0.7){
    hicrep_results %>%
    ggplot(
        aes(
            x=Genotype,
            y=hicrep.score,
            fill=is.downsampled,
            shape=ReadFilter
        )
    ) +
    geom_boxplot(outlier.size=0.5) +
    geom_jitter(
        data=. %>% filter(h.ideal == 'Ideal'),
        aes(
            x=Genotype,
            y=hicrep.score,
            fill=is.downsampled
        ),
        color='green',
        size=size
    ) +
    facet_grid2(
        rows=vars(Resolution),
        cols=vars(Celltype),
        scales='fixed'
    ) +
    theme(
        legend.position='top', 
        axis.text.x=element_text(angle=45, hjust=1)
    ) +
    add_ggtheme()
}
