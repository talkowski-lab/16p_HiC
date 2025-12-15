library(TADCompare)

###################################################
# Utilities
###################################################
TADCompare_load_matrix <- function(
    filepath,
    ...){
    load_mcool_file(
        filepath,
        type='df',
        cis=TRUE,
        ...
    ) %>% 
    select(c(range1, range2, IF))
}

load_hiTAD_results_for_TADCompare <- function(){
    check_cached_results(
        results_file=HITAD_TAD_RESULTS_FILE,
        force_redo=FALSE,
        results_fnc=load_all_hiTAD_TADs
    ) %>% 
    post_process_hiTAD_TAD_results() %>%
    filter(isMerged == 'Merged') %>% 
    mutate(resolution=scale_numbers(resolution, force_numeric=TRUE)) %>% 
    select(SampleID, resolution, chr, TAD.start, TAD.end) %>% 
    rename_with(~ str_replace(.x, 'TAD.', '')) %>% 
    mutate(region=chr) %>% 
    nest(TADs=c(chr, start, end)) %>% 
    rename('chr'=region)
}

load_cooltools_results_for_TADCompare <- function(){
    # Load boundary annotations
    check_cached_results(
        results_file=COOLTOOLS_TAD_RESULTS_FILE,
        force_redo=FALSE,
        results_fnc=load_all_cooltools_results,
        boundaries_only=TRUE
    ) %>% 
    # clean up 
    post_process_cooltools_results() %>% 
    mutate(resolution=scale_numbers(resolution, force_numeric=TRUE)) %>% 
    rename('TAD.params'=cooltools.params) %>% 
    select(resolution, TAD.params, SampleID, chr, bin.start) %>% 
    # remove entries with < 2 boundaries
    nest(boundaries=c(bin.start)) %>% 
    rowwise() %>% filter(nrow(boundaries) > 1) %>% 
    # convert boundaries to start/end format
    mutate(TADs=list(convert_boundaries_to_TADs(boundaries=boundaries))) %>% 
    ungroup() %>% select(-c(boundaries)) %>% unnest(TADs) %>% 
    # Nest for downstream analysis
    mutate(region=chr) %>% 
    group_by(resolution, TAD.params, SampleID, region) %>% 
    nest(TADs=c(chr, start, end)) %>% 
    ungroup() %>% 
    rename('chr'=region) %>% 
    select(resolution, TAD.params, SampleID, chr, TADs)
}

load_ConsensusTAD_results_for_TADCompare <- function(){
    check_cached_results(
        results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_ConsensusTAD_TADs
    ) %>% 
    # Only keep boundaries
    filter(isConsensusBoundary) %>% 
    # Clean up 
    post_process_ConsensusTAD_TAD_results() %>% 
    mutate(SampleID=glue('{Sample.Group}.Merged.Merged')) %>% 
    mutate(resolution=scale_numbers(resolution, force_numeric=TRUE)) %>% 
    select(resolution, TAD.params, SampleID, chr, bin.start) %>% 
    # remove entries with < 2 boundaries
    nest(boundaries=c(bin.start)) %>% 
    rowwise() %>% filter(nrow(boundaries) > 1) %>% 
    # convert boundaries to start/end format
    mutate(TADs=list(convert_boundaries_to_TADs(boundaries=boundaries))) %>% 
    ungroup() %>% select(-c(boundaries)) %>% unnest(TADs) %>% 
    # Nest for downstream analysis
    mutate(region=chr) %>% 
    group_by(resolution, TAD.params, SampleID, region) %>% 
    nest(TADs=c(chr, start, end)) %>% 
    ungroup() %>% 
    rename('chr'=region) %>% 
    select(resolution, TAD.params, SampleID, chr, TADs)
}

load_all_TAD_results_for_TADCompare <- function(){
    # hiTAD TAD results
    hiTAD.TADs.df <- 
        load_hiTAD_results_for_TADCompare() %>% 
        add_column(
            TAD.params=NULL,
            TAD.method='hiTAD'
        )
    # cooltools boundary results
    cooltools.TADs.df <- 
        load_cooltools_results_for_TADCompare() %>% 
        add_column(TAD.method='cooltools')
    # ConsensusTAD TAD results 
    consensusTAD.TADs.df <- 
        load_ConsensusTAD_results_for_TADCompare() %>% 
        add_column(TAD.method='ConsensusTAD')
    # default TADCompare method estimates TADs itself, include nothing
    spectralTAD.TADs.df <- 
        expand_grid(
            SampleID=unique(hiTAD.TADs.df$SampleID),
            chr=CHROMOSOMES,
            resolution=parsed.args$resolutions
        ) %>% 
        add_column(
            TADs=NULL,
            TAD.params=NULL,
            TAD.method='spectralTAD'
        )
    # Bind everything together
    bind_rows(
        hiTAD.TADs.df,
        cooltools.TADs.df,
        consensusTAD.TADs.df,
        spectralTAD.TADs.df
    ) %>%
    unite(
        'TAD.set.index',
        sep='~',
        remove=FALSE,
        c(
          TAD.method,
          TAD.params,
          resolution
        )
    ) %>% 
    rename('pre_tads'=TADs)
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

load_ConsensusTADs <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
    )
}

load_all_ConsensusTAD_TADs <- function(...){
    file.path(TAD_DIR, 'method_ConsensusTAD') %>%
    parse_results_filelist(
        suffix='-ConsensusTADs.tsv',
        filename.column.name='Sample.Group',
    ) %>% 
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
    rename('chr'=region)
}

###################################################
# Generate TADCompare results
###################################################
run_TADCompare <- function(
    filepath.Numerator,
    filepath.Denominator,
    Sample.Group.Numerator,
    Sample.Group.Denominator,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    pre_tads,
    ...){
    # row_index=1; filepath.Numerator=tmp$filepath.Numerator[[row_index]]; filepath.Denominator=tmp$filepath.Denominator[[row_index]]; Sample.Group.Numerator=tmp$Sample.Group.Numerator[[row_index]]; Sample.Group.Denominator=tmp$Sample.Group.Denominator[[row_index]]; resolution=tmp$resolution[[row_index]]; normalization=tmp$normalization[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]; z_thresh=tmp$z_thresh[[row_index]]; window_size=tmp$window_size[[row_index]]; gap_thresh=tmp$gap_thresh[[row_index]]; pre_tads=tmp$pre_tads[[row_index]];
    # chr1 @ 10Kb -> 24896x24896 matrix -> 40Gb is enough
    # Run TADCompare on the 2 matrices being compared
    matrix.numerator <-
        TADCompare_load_matrix(
            filepath.Numerator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    matrix.denominator <-
        TADCompare_load_matrix(
            filepath.Denominator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    tad.compare.results <- 
        TADCompare(
            matrix.numerator,
            matrix.denominator,
            resolution=resolution,
            z_thresh=z_thresh,
            window_size=window_size,
            gap_thresh=gap_thresh,
            pre_tads=pre_tads
        )
    # Format results to include boundary+gap scores for all bins + differential annotations
    tad.compare.results$Boundary_Scores %>% 
    as_tibble() %>%
    left_join(
        tad.compare.results$TAD_Frame %>%
        as_tibble() %>% 
        add_column(isTADBoundary=TRUE),
        suffix=c('.All', '.TADs'),
        by=join_by(Boundary)
    ) %>% 
    mutate(
        isTADBoundary=ifelse(is.na(isTADBoundary), FALSE, isTADBoundary),
        Enriched.Condition=
            case_when(
                Enriched_In.TADs == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.TADs == 'Matrix 2' ~ Sample.Group.Denominator,
                Enriched_In.All  == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.All  == 'Matrix 2' ~ Sample.Group.Denominator,
                TRUE                      ~ NA
            ),
        Differential=
            case_when(
                is.na(Differential.TADs) ~ Differential.All,
                TRUE                     ~ Differential.TADs
            ),
        Differential=
            ifelse(
                Differential %in% c('Differential', 'Non-Differential'),
                Differential, 
                'Differential'
            ),
        Type=
            case_when(
                is.na(Type.TADs) ~ Type.All,
                TRUE             ~ Type.TADs
            ),
        TAD_Score1=
            case_when(
                is.na(TAD_Score1.TADs) ~ TAD_Score1.All,
                TRUE                   ~ TAD_Score1.TADs
            ),
        TAD_Score2=
            case_when(
                is.na(TAD_Score2.TADs) ~ TAD_Score2.All,
                TRUE                   ~ TAD_Score2.TADs
            ),
        Gap_Score=
            case_when(
                is.na(Gap_Score.TADs) ~ Gap_Score.All,
                TRUE                  ~ Gap_Score.TADs
            )
    ) %>%
    rename(
        'TAD.Score.Numerator'=TAD_Score1,
        'TAD.Score.Denominator'=TAD_Score2,
        'TAD.isDifferential'=Differential,
        'TAD.Difference.Type'=Type
    ) %>% 
    rename_with(~ str_replace_all(.x, '_', '.')) %>% 
    select(-c(ends_with('.All'), ends_with('.TADs')))
}

run_all_TADCompare <- function(
    comparisons.df,
    hyper.params.df,
    TADs.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    ...){
    # chromosomes=CHROMOSOMES; force_redo=TRUE;
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    left_join(
        TADs.df,
        by=join_by(isMerged, resolution, chr)
    ) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TAD_DIR,
                'results_TADCompare',
                glue('merged_{isMerged}'),
                glue('method_{TAD.method}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('resolution_{scale_numbers(resolution)}'),
                # glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    arrange(resolution) %>% 
        {.}
    future_pmap(
    # pmap(
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

load_TADCompare_results <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
    )
}

load_all_TADCompare_results <- function(...){
    file.path(TAD_DIR, 'results_TADCompare') %>%
    parse_results_filelist() %>% 
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_TADCompare_results,
                .progress=TRUE
            )
    ) %>%
    unnest(TADs) %>% 
    select(-c(filepath))
}

post_process_TADCompare_results <- function(results.df){
    results.df %>%
    standardize_data_cols()
}

