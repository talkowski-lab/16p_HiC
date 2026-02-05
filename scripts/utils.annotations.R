###################################################
# Dependencies
###################################################
# needed for CTCF stuff
library(AnnotationHub)
library(GenomicRanges)
library(plyranges)

###################################################
# Load Specific Data
###################################################
load_sample_metadata <- function(filter=TRUE){
    SAMPLE_METADATA_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    {
        if(filter) {
            filter(., Included)
        } else {
            .
        }
    }
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=c('chr', 'chr.size.bp')
    )
}

get_min_resolution_per_matrix <- function(
    df=NULL,
    as_int=TRUE,
    filter_res=TRUE){
    # df=contacts.df; as_int=TRUE; filter_res=TRUE
    # get minimum viable resolution for each matrix based on Rao et at. 2014 definition
    MIN_SAMPLE_RESOLUTION_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    dplyr::select(SampleID, resolution) %>% 
    mutate(resolution=scale_numbers(resolution, force_numeric=as_int)) %>% 
    add_column(is.smallest.resolution=TRUE) %>% 
    {
        if (is.null(df)) {
            .
        } else {
            left_join(
                .,
                df,
                by=join_by(SampleID)
            )
        }
    } %>% 
    { 
        if (filter_res) {
            filter(., is.smallest.resolution) %>%
            dplyr::select(-c(is.smallest.resolution))
        } else {
            .
        }
    }
}

###################################################
# Gene Expression Data
###################################################
load_gene_expression_data <- function(){
    EXPRESSION_DATA_DIR %>%
    parse_results_filelist(
        suffix='-expression.tsv',
        filename.column.name='SampleID',
        pattern=NA
    ) %>%
    mutate(
        expr=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                show_col_types=FALSE,
                progress=FALSE,
                .progress=TRUE
            )
    ) %>%
    unnest(expr) %>% 
    mutate(chr=rename_chrs(chr, to_label=TRUE)) %>% 
    select(
        -c(
            starts_with('X16p'),
            starts_with('GM_'),
            geneid,
            filepath
        )
    ) %>% 
    rename(
        'EnsemblID'=ensemblid,
        'TPM.mean'=mean_tpm,
        'TPM.sd'=sd_tpm
    )
}

###################################################
# Genome Feature Annotations
###################################################
load_CTCF_sites <- function(force.redo=FALSE){
    # https://dozmorovlab.github.io/CTCF/
    # annotation object
    check_cached_results(
        results_file=CTCF_SITE_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                AnnotationHub() %>% 
                # subset specific set of CTCF sites predicted in humans by JASPAR 2022
                subset(
                    species == 'Homo sapiens' &
                    genome == 'hg38' &
                    dataprovider == 'JASPAR 2022' &
                    preparerclass == "CTCF"
                ) %>%
                {.[['AH104727']]} %>%
                # make tidy tibble
                as_tibble() %>% 
                # separate_wider_delim(
                #     name,
                #     delim='_',
                #     names=
                #         c(
                #             'db',
                #             'db.set',
                #             'bioset',
                #             'genomic.feature',
                #             'MotifID'
                #         )
                # ) %>%
                dplyr::rename(
                    'length'=width,
                    'chr'=seqnames
                ) %>% 
                select(
                    -c(
                        # db,
                        # db.set,
                        # bioset,
                        # genomic.feature,
                        sequence
                    )
                )
            }
    )
}

