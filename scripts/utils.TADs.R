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
    ) %>% 
    dplyr::rename('bin'=start) %>% 
    select(-c(end))
}

load_DI_data <- function(
    filepath,
    method,
    ...){
    if (method == 'hiTAD') {
        load_hiTAD_DIs(filepath, ...)
    } else if (method == 'cooltools') {
        load_cooltools_DIs(filepath, ...)
    } else {
        stop(glue('Invalid method ({method})'))
    }
}

annotated_bins_with_TADs <- function(
    DIs,
    TADs,
    ...){
    # DIs=tmp$DIs[[1]]; TADs=tmp$TADs[[1]]
    # DIs %>% nrow(); TADs %>% nrow()
    DIs %>% 
    # Annotated within TAD bins
    left_join(
        TADs %>% 
        select(TAD.start, TAD.end) %>% 
        add_column(is.TAD.interior=TRUE),
        by=join_by(between(bin, TAD.start, TAD.end, bounds="()"))
    ) %>%
    # Annotated within TAD Start 
    left_join(
        TADs %>% 
        select(TAD.start) %>% 
        add_column(is.TAD.start=TRUE),
        by=join_by(bin == TAD.start)
    ) %>%
    # Annotated within TAD Ends
    left_join(
        TADs %>% 
        select(TAD.end) %>% 
        add_column(is.TAD.end=TRUE),
        by=join_by(bin == TAD.end)
        # by=join_by(between(bin, TAD.start, TAD.end))
    ) %>% 
    mutate(across(starts_with('is.TAD'), ~ !is.na(.x))) %>% 
    mutate(
        TAD.status=
            case_when(
                 is.TAD.start &  is.TAD.end ~ 'TAD Split',
                 is.TAD.start & !is.TAD.end ~ 'TAD Start',
                !is.TAD.start &  is.TAD.end ~ 'TAD End',
                 is.TAD.interior            ~ 'Inside TAD',
                 TRUE                       ~ 'Outside TADs'
            )
    ) %>%
    select(bin, end, DI, TAD.status)
}

list_all_DI_results <- function(
    pattern=NULL,
    ignore.case=TRUE,
    invert=FALSE){

    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-DI.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(weight == 'ICE') %>%
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
    standardize_data_cols()
}

load_all_DI_data <- function(
    hitad.TAD.df,
    ...){
    list_all_DI_results(...) %>% 
    mutate(
        DIs=
            pmap(
                .,
                load_DI_data,
                .progress=TRUE
            )
    ) %>%
    unnest(DIs) %>% 
    standardize_data_cols() %>% 
    select(-c(filepath))
}

annotate_DI_with_TADs <- function(
    hitad.DI.df,
    hitad.TAD.df,
    pair_columns=
        c(
            'isMerged',
            'resolution',
            'SampleID',
            'chr'
        ),
    ...){
    hitad.DI.df %>% 
    nest(DIs=c(bin, DI)) %>% 
    post_process_hiTAD_DI_results() %>% 
    # Match TAD annotations to bin-wise DI results
    inner_join(
        hitad.TAD.df %>% 
        nest(TADs=c(TAD.start, TAD.end, TAD.length)),
        by=pair_columns,
    ) %>%
    # Annotate every bin as to whether it is a TAD boundary/inside/outside
    mutate(
       tad.annotated.bins=
            pmap(
                .l=.,
                .f=annotated_bins_with_TADs,
                .progress=TRUE
            )
    ) %>% 
    unnest(tad.annotated.bins) %>%
    select(-c(TADs, DIs))
}

list_all_dataset_pairs <- function(
    pair_grouping_cols,
    ...){
    # pair_grouping_cols=c('isMerged', 'resolution')
    list_all_DI_results(...) %>% 
    standardize_data_cols() %>% 
    nest(
        SampleInfo=
            c(
                # weight,
                SampleID,
                Edit,
                Genotype,
                Celltype,
                CloneID,
                TechRepID
            )
    ) %>% 
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
    select(
        all_of(
            c(
                pair_grouping_cols,
                'SamplePairInfo',
                'filepath.P1',
                'filepath.P2'
            )
        )
    ) %>% 
    unnest(SamplePairInfo)
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

post_process_hiTAD_TAD_results <- function(results.df){
    results.df %>% 
    # Subset to only relevant parameters
    filter(weight == 'ICE') %>% 
    standardize_data_cols() %>% 
    mutate(TAD.length=TAD.end - TAD.start) %>% 
    select(
        -c(
            method,
            weight,
            threshold,
            mfvp,
            ReadFilter
        )
    )
}

post_process_cooltools_TAD_results <- function(results.df){
    results.df %>% 
    standardize_data_cols() %>% 
    mutate(
        window.size.bins=window.size / resolution,
        window.size=scale_numbers(window.size),
        param.combo=glue('{weight},{mfvp},{window.size.bins},{threshold}')
    ) %>% 
    rename('TAD.boundary'=TAD.start) %>% 
    select(-c(TAD.end, method))
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

calculate_profile_differences <- function(
    joint.DIs,
    group_cols=NULL){
    # Compare DI profile of 2 samples and calc correlation/diff
    joint.DIs %>% 
    {
        if (is.null(group_cols)) {
            .
        } else {
            group_by(., across(all_of(group_cols))) %>%
            dplyr::add_count() %>%
            filter(n > 2) %>%
            select(-c(n))
        }
    } %>% 
    summarize(
        pearson=
            cor.test(
                DI.P1, DI.P2,
                alternative='two.sided',
                exact=FALSE,  # avoid "ties" error
                method='pearson'
            ) %>% 
            broom::tidy() %>%
            rename('corr'=estimate),
        spearman=
            cor.test(
                DI.P1, DI.P2,
                alternative='two.sided',
                exact=FALSE,  # avoid "ties" error
                method='spearman'
            ) %>% 
            broom::tidy() %>% 
            rename('corr'=estimate),
        paired_t=
            t.test(
                DI.P1, DI.P2,
                paired=TRUE,
                exact=FALSE,
                mu=0,
                var.equal=FALSE,
                alternative="two.sided"
            ) %>% 
            broom::tidy(),
        paired_wilcox=
            wilcox.test(
                DI.P1, DI.P2,
                paired=TRUE,
                exact=FALSE,
                alternative="two.sided"
            ) %>% 
            broom::tidy()
    ) %>%
    ungroup() %>% 
    {
        if (is.null(group_cols)) {
            pivot_longer(
                .,
                everything(),
                names_to='tmp',
                values_to='value'
            )
        } else {
            pivot_longer(
                .,
                -c(all_of(group_cols)),
                names_to='tmp',
                values_to='value'
            )
        }
    } %>% 
    unnest(c(value)) %>% 
    {
        if (is.null(group_cols)) {
            select(., tmp, corr, statistic, p.value)
        } else {
            select(., all_of(group_cols), tmp, corr, statistic, p.value) %>% 
            ungroup()
        }
    } %>% 
    dplyr::rename('method'=tmp)
}

calculate_pair_InsulationCorr <- function(
    DIs.P1,
    DIs.P2,
    # group_cols=c('TAD.status'),
    ...){
    # DIs.P1=tmp$DIs.P1[[1]]; DIs.P2=tmp$DIs.P2[[1]]; # group_cols=c('TAD.status')
    # Join tracks together by bin
    joint.DIs <- 
        inner_join(
            DIs.P1,
            DIs.P2,
            suffix=c('.P1', '.P2'),
            by=join_by(bin)
        )
    # Compute correlation stats over all bins
    total.results <- 
        joint.DIs %>%
        calculate_profile_differences()
    # Compute correlation stats over all bins
    group.results <- 
        joint.DIs %>%
        filter(TAD.status.P1 == TAD.status.P2) %>% 
        calculate_profile_differences(
            group_cols=c('TAD.status.P1', 'TAD.status.P2')
        )
    # Bind results
    bind_rows(
        total.results %>% 
            add_column(TAD.statuses='all'),
        group.results %>%
            unite(
                'TAD.statuses',
                # all_of(group_cols),
                c('TAD.status.P1', 'TAD.status.P2'),
                sep="#",
                remove=FALSE
            )
    )
}

calculate_all_pairs_InsulationCorr <- function(
    DI.annotations,
    pair_grouping_cols,
    ...){
    # tad.annotations Must contain the following columns
    # SampleInfo: nested tibble of Sample attributes
    # DIs: nested tibble with 2+ columns; bin, DI, bin_info...
    # pair_grouping_cols: all columns listed here, only compare pairs of TAD sets that match
    # that match along these attributes e.g. resolution, chr
    # pair_grouping_cols=c('isMerged', 'resolution', 'chr'); DI.annotations=hitad.DI.df %>% nest(DIs=c(bin, DI, TAD.status)) %>% nest(SampleInfo=c(Edit, Genotype, Celltype, CloneID, TechRepID, SampleID)); DI.annotations
    DI.annotations %>% 
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
        results=
            pmap(
                .l=.,
                .f=calculate_pair_InsulationCorr,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>%
    select(-c(TAD.status.P1, DIs.P1, TAD.status.P2, DIs.P2))
}

###################################################
# Plot stuff
###################################################
