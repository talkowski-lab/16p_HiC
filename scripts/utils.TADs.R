library(tidyverse)
library(stringi)
# library(HiContacts)
library(glue)
library(furrr)
library(ggplot2)
library(viridis)
library(ggh4x)
# library(ggtext)
library(cowplot)
###################################################
# Misc
###################################################
load_cooltools_file <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        progress=FALSE,
        show_col_types=FALSE
    ) %>%
    rename(
        'chr'=chrom,
        'TAD.start'=start,
        'TAD.end'=end
    ) %>% 
    pivot_longer(
        ends_with('000'),
        names_to='TAD.stat',
        values_to='value'
    ) %>%
    # Need to reverse since some stat names have _ in them
    mutate(TAD.stat=stri_reverse(TAD.stat)) %>% 
    separate_wider_delim(
        TAD.stat,
        delim='_',
        names=c(
            'window.size',
            'stat'
        ),
        too_many='merge'
    ) %>%
    mutate(
        across(
            c(stat, window.size),
            stri_reverse
        )
    ) %>% 
    pivot_wider(
        names_from=stat,
        values_from=value
    ) %>%
    mutate(
        is_boundary=as.logical(is_boundary),
        window.size=as.integer(window.size)
    )
}

###################################################
# Load Insulation data
###################################################
load_cooltools_DIs <- function(
    filepath,
    ...){
    filepath %>% 
    load_cooltools_file() %>% 
    # filter(!is_bad_bin) %>% 
    select(
        is_bad_bin,
        # is_boundary,
        log2_insulation_score,
        # n_valid_pixels,
        # boundary_strength,
        # sum_balanced,
        # sum_counts,
        region,
        chrom,
        start,
        end
    )
}

load_hiTAD_DIs <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'start',
                'end',
                'DI'
            )
    )
}

load_DI_data <- function(
    filepath,
    method,
    filetype,
    ...){
    if (method == 'hiTAD' & filetype == 'DI') {
        load_hiTAD_DIs(filepath, ...)
    } else if (method == 'cooltools' & filetype == 'TAD') {
        load_cooltools_DIs(filepath, ...)
    } else {
        message(glue('Invalid method ({method}) and/or suffix ({suffix})'))
    }
}

load_all_DI_data <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-(TADs|DI).tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    process_matrix_name() %>% 
    mutate(
        filetype=
            filepath %>%
            str_extract('-(TAD|DI).tsv$', group=1),
        DIs=
            pmap(
                .,
                load_DI_data,
                .progress=TRUE
            )
    ) %>%
    select(-c(filepath)) %>% 
    unnest(DIs)
}

###################################################
# Load TAD annotation data
###################################################
load_arrowhead_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        skip=2,
        show_col_types=FALSE,
        col_select=c(1,2,3,12,13,14),
        # col_select=c(1,2,3,4,5,6,12,13,14),
        col_names=
            c(
                'chr',
                'start',
                'end',
                'B.chr',
                'B.start',
                'B.end',
                'name',
                'score',
                'strand1',
                'strand2',
                'color',
                'corner.score',
                'uVarScore',
                'lVarScore',
                'upSign',
                'loSign'
            )
    )
}

load_hiTAD_TAD_annotation <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
        col_names=
            c(
                'chr',
                'TAD.start',
                'TAD.end'
            )
    )
}

load_cooltools_TAD_annotation <- function(
    filepath,
    ...){
    # filepath="/data/talkowski/Samples/16p_HiC/results/TADs/method_cooltools/resolution_100000/weight_ICE/threshold_0/mfvp_0.33/16p.DEL.A3.NSC.HiC.hg38.mapq_30.1000-TAD.tsv"
    filepath %>% 
    load_cooltools_file() %>% 
    filter(!is_bad_bin) %>% 
    filter(is_boundary) %>% 
    select(
        # is_bad_bin,
        # is_boundary,
        # log2_insulation_score,
        n_valid_pixels,
        boundary_strength,
        window.size,
        # sum_balanced,
        # sum_counts,
        region,
        chr,
        TAD.start,
        TAD.end
    )
}

load_TAD_annotation <- function(
    filepath,
    method,
    ...){
    if (method == 'hiTAD') {
        load_hiTAD_TAD_annotation(filepath)
    } else if (method == 'cooltools') {
        load_cooltools_TAD_annotation(filepath)
    } else if (method == 'arrowhead') {
        load_arrowhead_TAD_annotation(filepath)
    } else {
        message(glue('Invalid method {method}'))
    }
}

load_all_TAD_annotations <- function(
    pattern=NULL,
    ignore.case=TRUE,
    invert=FALSE){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-TAD.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    {
        if (!is.null(pattern)) {
            if (invert) {
                filter(., !grepl(pattern, filepath, ignore.case=ignore.case))
            } else {
                filter(.,  grepl(pattern, filepath, ignore.case=ignore.case))
            }
        } else {
            .
        }
    } %>% 
    mutate(
        TADs=
            # future_pmap(
            pmap(
                .,
                load_TAD_annotation,
                .progress=TRUE
            )
    ) %>%
    unnest(TADs) %>% 
    select(-c(filepath))
}

###################################################
# Compute Similarities
###################################################
define_TAD_pairs <- function(
    TADs.P1,
    TADs.P2,
    ...){
    # TADs.P1=tmp$TADs.P1[[1]]; TADs.P2=tmp$TADs.P2[[1]];
    cross_join(
        TADs.P1 %>% 
            unite(
                'TAD_idx',
                TAD.start, TAD.end,
                sep='#',
                remove=FALSE
            ),
        TADs.P2 %>% 
            unite(
                'TAD_idx',
                TAD.start, TAD.end,
                sep='#',
                remove=FALSE
            ),
        suffix=c('.P1', '.P2')
    ) %>%
    # only keep pairs of TADs that overlap at all
    filter(
        between(TAD.start.P1, TAD.start.P2, TAD.end.P2) |
        between(TAD.end.P1,   TAD.start.P2, TAD.end.P2) |
        between(TAD.start.P2, TAD.start.P1, TAD.end.P1) |
        between(TAD.end.P2,   TAD.start.P1, TAD.end.P1)
    ) %>% 
    #   000000000111111111122222222223333
    #   123456789012345678901234567890123
    # P ---|++++|------|+++++|-----------|+++++++++++|------|++++++|------------|+++++|--
    # Q ------|++++|------------|++++|-------|++++|------|+++++++++++++++|----|++++|-----
    # P 04-09, 16-22
    # Q 07-12, 25-30
    mutate(
        rightmost.start=max(TAD.start.P1, TAD.start.P2),
        leftmost.end=min(TAD.end.P1, TAD.end.P2),
        overlap=leftmost.end  - rightmost.start,
        # observed overlap / total possible overlap
        TAD.jaccard=((overlap**2) / (TAD.length.P1 * TAD.length.P2))
    ) %>%
    select(-c(rightmost.start, leftmost.end, overlap))
}

calculate_pair_MoC <- function(
    TADs.P1,
    TADs.P2,
    ...){
    # TADs.P1=tmp$TADs.P1[[1]]; TADs.P2=tmp$TADs.P2[[1]]
    # TADs.P1, TADs.P2 should be tibble with columns TAD.start, TAD.end, TAD.length
    # 1 row per TAD
    # Number of TADs
    nTADs.P <- nrow(TADs.P1)
    nTADs.Q <- nrow(TADs.P2)
    # Skip trivial cases
    if (nTADs.P == nTADs.Q & nTADs.P == 1) {
        1
    } else {
        norm_const <- 1 / (sqrt(nTADs.P * nTADs.Q) - 1)
        define_TAD_pairs(
            TADs.P1,
            TADs.P2
        ) %>% 
        pull(TAD.jaccard) %>% 
        sum() %>% 
        {(. - 1) * norm_const}
    }
}

calculate_all_pairs_MoCs <- function(
    tad.annotations,
    pair_grouping_cols,
    ...){
    # pair_grouping_cols=c('isMerged', 'resolution', 'chr'); tad.annotations=hitad.annotations %>% nest(TADs=c(TAD.start, TAD.end, TAD.length)) %>% nest(SampleInfo=c(Edit, Genotype, Celltype, CloneID, TechRepID, weight, SampleID)); tad.annotations
    # tad.annotations Must contain the following columns
    # SampleInfo: nested tibble of Sample attributes
    # TADs: nested tibble with 3 columns; TAD.start, TAD.end, TAD.length
    # pair_grouping_cols: all columns listed here, only compare pairs of TAD sets that match
    # that match along these attributes e.g. resolution, chr
    tad.annotations %>% 
    # List all pairs of annotations (SampleID + chr) that also have matching parameter values listed
    get_all_row_combinations(
        cols_to_pair=pair_grouping_cols,
        keep_self=FALSE
    ) %>% 
    rowwise() %>% 
    mutate(
        SamplePairInfo=
            merge_sample_info(
                SampleInfo.P1,
                SampleInfo.P2
            ) %>%
            list()
    ) %>% 
    ungroup() %>% 
    select(-c(starts_with('SampleInfo'))) %>% 
    unnest(SamplePairInfo) %>% 
    # Finally compute all MoCs for all listed pairs of annotations
    mutate(
        mocs=
            pmap(
                .l=.,
                .f=calculate_pair_MoC,
                .progress=TRUE
            )
    ) %>%
    unnest(mocs) %>%
    select(-c(TADs.P1, TADs.P2))
}

calculate_pair_boundary_conservation <- function(
    TADs.P1,
    TADs.P2,
    tolerance=50000,
    ...){
    # TADs.P1, TADs.P2 should be tibble with columns TAD.start, TAD.end, TAD.length, resolution
    # 1 row per TAD
    cross_join(
        boundaries.P1,
        boundaries.P2,
        suffix=c('.P', '.Q')
    )
}

###################################################
# Plot stuff
###################################################
