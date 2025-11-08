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
# Cooltools 
###################################################
load_cooltools_file <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        progress=FALSE,
        show_col_types=FALSE
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
    )
}

load_all_cooltools_results <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-TAD.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(method == 'cooltools') %>% 
    mutate(
        DIs=
            # pmap(
            future_pmap(
                .,
                load_cooltools_file,
                .progress=TRUE
            )
    ) %>%
    unnest(DIs) %>% 
    select(-c(filepath))
}

post_process_cooltools_results <- function(results.df){
    results.df %>% 
    filter(!is_bad_bin) %>% 
    standardize_data_cols() %>% 
    rename(
        'chr'=chrom,
        'bin.start'=start
    ) %>% 
    mutate(
        window.size.bins=window.size / resolution,
        window.size=scale_numbers(window.size),
        is.boundary=as.logical(is.boundary)
    ) %>% 
    unite(
        'cooltools.params',
        weight, window.size.bins, mfvp, threshold, 
        sep='#',
        remove=TRUE
    ) %>% 
    select(
        -c(
            method,
            end,
            is_bad_bin,
            sum_counts,           
            sum_balanced,         
            ReadFilter,
            region,               
            n_valid_pixels       
            # method               
            # resolution           
            # weight               
            # threshold            
            # mfvp                 
            # window.size          
            # Edit                 
            # Celltype             
            # Genotype             
            # CloneID              
            # TechRepID            
            # isMerged             
            # SampleID             
            # chr                
            # bin.start                
            # is.bad.bin
            # log2.insulation.score
            # boundary.strength    
            # is.boundary          
        )
    ) %>% 
    rename_with(
        c(
            log2_insulation_score,
            boundary_strength,
            is_boundary
        ),
        ~ str_replace_all(.x, '_', '.')
    )
}

###################################################
# hiTAD 
###################################################
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

load_all_hiTAD_DIs <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-DI.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(method == 'hiTAD') %>% 
    filter(weight == 'ICE') %>%
    mutate(
        DIs=
            # pmap(
            future_pmap(
                .,
                load_hiTAD_DIs,
                .progress=TRUE
            )
    ) %>%
    unnest(DIs) %>% 
    select(-c(filepath))
}

load_hiTAD_TADs <- function(
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

load_all_hiTAD_TADs <- function(){
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-TAD.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(method == 'hiTAD') %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_hiTAD_TADs,
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

###################################################
# Data loading wrappers
###################################################
list_all_dataset_pairs <- function(
    pair_grouping_cols,
    ...){
    # pair_grouping_cols=c('isMerged', 'resolution')
    parse_results_filelist(
        input_dir=TAD_DIR,
        suffix='-DI.tsv'
    ) %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(
        (method == 'hiTAD' & weight == 'ICE') |
        (method == 'cooltools')
    ) %>%
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
        cols_to_pair=c(pair_grouping_cols, 'method'),
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
# Compute Similarities
###################################################
define_TAD_pairs <- function(
    TADs.P1,
    TADs.P2,
    ...){
    # TADs.P1=tmp$TADs.P1[[1]]; TADs.P2=tmp$TADs.P2[[1]];
    cross_join(
        TADs.P1,
        TADs.P2,
        suffix=c('.P1', '.P2')
    ) %>%
    # only keep pairs of TADs that overlap at all
    #   000000000111111111122222222223333
    #   123456789012345678901234567890123
    # P ---|++++|------|+++++|-----------|+++++++++++|------|++++++|------------|+++++|--
    # Q ------|++++|------------|++++|-------|++++|------|+++++++++++++++|----|++++|-----
    filter(
        between(TAD.start.P1, TAD.start.P2, TAD.end.P2) |
        between(TAD.end.P1,   TAD.start.P2, TAD.end.P2) |
        between(TAD.start.P2, TAD.start.P1, TAD.end.P1) |
        between(TAD.end.P2,   TAD.start.P1, TAD.end.P1)
    ) %>% 
    # Only keep each TAD once at most
    # compute overlap similarity (MoC) for each pair of TADs between the 2 annotations
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
    # pair_grouping_cols=c('isMerged', 'resolution', 'chr'); tad.annotations=hitad.TAD.df %>% nest(TADs=c(TAD.start, TAD.end, TAD.length)) %>% nest(SampleInfo=c(Edit, Genotype, Celltype, CloneID, TechRepID, SampleID)); tad.annotations
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
        # {.} -> tmp
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
# TADCompare
###################################################
set_up_sample_groups <- function(sample.groups){
    list_mcool_files() %>%
    get_min_resolution_per_matrix() %>% 
    distinct() %>% 
    # Now group samples by condition, 
    filter(!isMerged) %>% 
    nest(samples.df=-c(isMerged)) %>% 
    cross_join(sample.groups) %>%
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.Pattern) ~ Sample.Group,
                        TRUE ~ NA
                    )
            ) %>%
            filter(!is.na(Sample.Group)) %>% 
            list()
    ) %>%
    # minimum and max resoltion of all individual matrices per comparison
    mutate(
        resolution.min=min(samples.df$resolution),
        resolution.max=max(samples.df$resolution)
    ) %>% 
    ungroup() %>%
    mutate(resolution=list(unique(c(resolution.min, resolution.max)))) %>%
    unnest(resolution) %>% 
    mutate(
        resolution.type=
            case_when(
                resolution == resolution.max ~ 'max',
                resolution == resolution.min ~ 'min',
                TRUE                         ~ NA
            )
    ) %>% 
    select(
        -c(
            resolution.min,
            resolution.max,
            ends_with('.Pattern')
        )
    )
    
}

run_ConsensusTAD <- function(
    samples.df,
    resolution,
    normalization,
    range1,
    range2,
    ...){
    samples.df %>%
    pull(filepath) %>% 
    sapply(
        load_mcool_file,
        resolution=resolution,
        normalization=normalization,
        range1=range1,
        range2=range2,
        cis=TRUE,
        type='dense',
        simplify=FALSE
    ) %>% 
    ConsensusTADs(
        resolution=resolution,
        ...
    ) %>%
    {.$All_Regions}
}

run_all_ConsensusTADs <- function(
    sample.groups.df,
    hyper.params.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    ...){
    sample.groups.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TAD_DIR,
                'results_ConsensusTADs',
                glue('merged_{isMerged}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group}-ConsensusTADs.tsv')
            )
    ) %>% 
    arrange(resolution) %>% 
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_ConsensusTAD,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

run_TADCompare <- function(
    matrix1.filepath,
    matrix2.filepath,
    matrix1.name,
    matrix2.name,
    resolution,
    normalization,
    range1,
    range2,
    ...){
    # row_index=16; samples.df=tmp$samples.df[[row_index]]; sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL; resolution=tmp$resolution[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]; effect.col='Sample.Group'; p.method='fdr';
    matrix1 <-
        load_mcool_file(
            matrix1.filepath,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2,
            cis=TRUE,
            type='dense'
        )
    matrix2 <-
        load_mcool_file(
            matrix2.filepath,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2,
            cis=TRUE,
            type='dense'
        )
    TADCompare(
        matrix1,
        matrix2,
        resolution=resolution,
        ...
    ) %>%
    {.$Boudary_Scores} %>% 
    mutate(
        Enriched_Condition=
            case_when(
                Enriched_In == 'Matrix 1' ~ matrix1.name,
                Enriched_In == 'Matrix 2' ~ matrix2.name,
                TRUE                     ~ NA
            )
    )
}

run_all_TADCompare <- function(
    comparisons.df,
    hyper.params.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    sample_group_priority_fnc,
    ...){
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TAD_DIR,
                'results_TADCompare',
                glue('merged_{isMerged}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    arrange(resolution) %>% 
        {.} -> tmp
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_TADCompare,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

