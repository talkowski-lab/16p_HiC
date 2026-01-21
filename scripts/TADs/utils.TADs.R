library(tidyverse)
library(stringi)
library(glue)
library(furrr)
library(TADCompare)

###################################################
# Utilities
###################################################
convert_boundaries_to_TADs <- function(
    boundaries,
    start.col.name='start',
    end.col.name='end',
    ...){
    # convert list of TAD boundaries to start/end format
    # every bnoundary is a start+end excpet for the first and last ones
    boundaries <- unlist(boundaries)
    tibble(
        TAD.start=boundaries[1:length(boundaries)-1],
        TAD.end=boundaries[2:length(boundaries)]
    ) %>%
    add_row(TAD.start=NA, TAD.end=NA) %>% 
    rename(
       !!start.col.name := TAD.start,
       !!end.col.name := TAD.end
    )
}

convert_TADs_to_boundaries <- function(
    # convert list of TAD boundaries to start/end format
    # every bnoundary is a start+end excpet for the first and last ones
    TAD.starts,
    TAD.ends,
    boundary.col.name='boundary',
    ...){
    c(TAD.starts, TAD.ends) %>%
    unique() %>%
    sort() %>% 
    tibble(boundary=.) %>% 
    rename(boundary.col.name=boundary)
}

###################################################
# Cooltools 
###################################################
load_cooltools_file <- function(
    filepath,
    boundaries_only=FALSE,
    ...){
    # filepath=tmp$filepath[[1]]; boundaries_only=TRUE
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
    ) %>%
    rename('is.boundary'=is_boundary) %>% 
    mutate(is.boundary=as.logical(is.boundary)) %>% 
    {
        if (boundaries_only){
            filter(., is.boundary)
        } else {
            .
        }
    }
}

load_all_cooltools_results <- function(boundaries_only=FALSE){
    file.path(TAD_DIR) %>% 
    parse_results_filelist(suffix='-TAD.tsv') %>%
    filter(method == 'cooltools') %>% 
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    # {.} -> tmp
    mutate(
        DIs=
            # pmap(
            future_pmap(
                .,
                load_cooltools_file,
                boundaries_only=boundaries_only,
                .progress=TRUE
            )
    ) %>%
    unnest(DIs) %>% 
    select(-c(filepath))
}

post_process_cooltools_results <- function(results.df){
    results.df %>% 
    filter(!is_bad_bin) %>% 
    rename(
        'chr'=chrom,
        'bin.start'=start
    ) %>% 
    mutate(
        window.size.bins=window.size / resolution,
        window.size=scale_numbers(window.size, force_numeric=TRUE),
    ) %>% 
    unite(
        'cooltools.params',
        weight, window.size.bins, mfvp, threshold, 
        sep='#',
        remove=TRUE
    ) %>% 
    select(
        -c(
            # method,
            end,
            is_bad_bin,
            sum_counts,           
            sum_balanced,         
            n_valid_pixels,
            ReadFilter,
            region,               
            # Edit                 
            # Celltype             
            # Genotype             
            # CloneID              
            # TechRepID            
            # isMerged             
            # SampleID             
            # resolution           
            # weight               
            # threshold            
            # mfvp                 
            # window.size          
            # chr                
            # bin.start                
        )
    ) %>% 
    dplyr::rename_with(
        ~ str_replace_all(.x, '_', '.'),
        c(
            # is_bad_bin,
            # n_valid_pixels,
            # sum_counts,           
            # sum_balanced,         
            log2_insulation_score,
            boundary_strength
        )
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
                'start',
                'end'
            )
    )
}

load_all_hiTAD_TADs <- function(){
    TAD_DIR %>% 
    parse_results_filelist(suffix='-TAD.tsv') %>%
    filter(method == 'hiTAD') %>% 
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
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
    mutate(length=end - start) %>% 
    mutate(Sample.Group=str_replace_all(SampleID, '.Merged.Merged', '')) %>% 
    filter(weight == 'ICE') %>% 
    filter(isMerged) %>% 
    dplyr::select(
        -c(
            weight,
            threshold,
            mfvp,
            ReadFilter,
            isMerged,
            Edit, Celltype, Genotype, CloneID, TechRepID,
            SampleID
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
# Generate ConsensusTADs
###################################################
run_ConsensusTAD <- function(
    samples.df,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    ...){
    # row_index=750; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; normalization=tmp$normalization[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]
    # filepath=samples.df$filepath[[1]]
    sampleID.mapping <- 
        samples.df %>%
        mutate(
            SampleID=as.character(glue('{SampleID}.score')),
            og.sample=as.character(glue('Sample {row_number()}'))
        ) %>%
        select(SampleID, og.sample) %>% 
        deframe()
    # Generate ConsensusTAD results
    consensus.results <- 
        samples.df %>%
        pull(filepath) %>% 
        sapply(
            TADCompare_load_matrix,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2,
            simplify=FALSE,
            USE.NAMES=FALSE
        ) %>% 
        ConsensusTADs(
            resolution=resolution,
            z_thresh=z_thresh,
            window_size=window_size,
            gap_thresh=gap_thresh
        )
    # Save data for all bins + consensus boundaries
    consensus.results$All_Regions %>% 
    tibble() %>%
    # Specify which bins were actually called as boundaries
    left_join(
        consensus.results$Consensus %>% 
        tibble() %>% 
        select(Coordinate) %>% 
        add_column(isConsensusBoundary=TRUE),
        by=join_by(Coordinate)
    ) %>% 
    mutate(isConsensusBoundary=ifelse(is.na(isConsensusBoundary), FALSE, isConsensusBoundary)) %>% 
    rename(
        all_of(sampleID.mapping),
        'bin.start'=Coordinate,
        'consensus.score'=Consensus_Score
    )
}

run_all_ConsensusTADs <- function(
    sample.groups.df,
    hyper.params.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    ...){
    # chromosomes=CHROMOSOMES; force_redo=TRUE;
    sample.groups.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output files accordingly
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TAD_DIR,
                'method_ConsensusTAD',
                glue('merged_{isMerged}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group}-ConsensusTADs.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !(file.exists(results_file)))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            dplyr::count(
                .,
                z_thresh,
                window_size,
                gap_thresh,
                resolution,
                Sample.Group
            )
        )
    } %>%
    arrange(resolution) %>% 
    # future_pmap(
    pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_ConsensusTAD
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

load_ConsensusTADs <- function(filepath, ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
    )
}

load_all_ConsensusTAD_TADs <- function(...){
    TAD_DIR %>% 
    parse_results_filelist(
        suffix='-ConsensusTADs.tsv',
        filename.column.name='Sample.Group',
    ) %>% 
    filter(method == 'ConsensusTAD') %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_ConsensusTADs,
                .progress=TRUE
            )
    ) %>%
    unnest(TADs) %>% 
    select(-c(filepath))
}

post_process_ConsensusTAD_TAD_results <- function(results.df){
    results.df %>%
    unite(
        'TAD.params',
        sep='#',
        z.thresh,
        window.size,
        gap.thresh,
    ) %>% 
    # nest(scores=ends_with('.score'))
    # Only keep boundaries
    filter(isConsensusBoundary) %>% 
    # Clean up 
    rename('chr'=region) %>% 
    select(method, resolution, TAD.params, Sample.Group, chr, bin.start) %>% 
    nest(boundaries=c(bin.start)) %>% 
    # remove entries with < 2 boundaries
    rowwise() %>% filter(nrow(boundaries) > 1) %>% 
    # convert boundaries to start/end format
    mutate(
        TADs=
            list(
                convert_boundaries_to_TADs(
                    boundaries=boundaries,
                    start.col.name='start',
                    end.col.name='end',
                )
            )
    ) %>% 
    ungroup() %>% select(-c(boundaries)) %>% unnest(TADs) %>% 
    # Nest for downstream analysis
    mutate(length=end - start)
}

###################################################
# Compute Similarities
###################################################
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

