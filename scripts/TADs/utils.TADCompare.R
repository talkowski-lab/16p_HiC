library(TADCompare)

###################################################
# Load TADs to compare
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
        # force_redo=TRUE,
        results_fnc=load_all_hiTAD_TADs
    ) %>% 
    filter(weight == 'balanced') %>% 
    filter(isMerged) %>% 
    select(SampleID, resolution, chr, start, end) %>% 
    mutate(Sample.Group=str_replace_all(SampleID, '.Merged.Merged', '')) %>% 
    mutate(region=chr) %>% 
    mutate(length=end - start) %>% 
    nest(TADs=c(chr, start, end, length)) %>% 
    dplyr::rename('chr'=region)
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
    dplyr::rename('chr'=region) %>% 
    select(resolution, TAD.params, SampleID, chr, TADs)
}

load_ConsensusTAD_results_for_TADCompare <- function(){
    check_cached_results(
        results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_ConsensusTAD_TADs
    ) %>% 
    post_process_ConsensusTAD_TAD_results() %>% 
    mutate(SampleID=glue('{Sample.Group}.Merged.Merged')) %>% 
    dplyr::select(
        method, TAD.params, resolution, 
        SampleID, Sample.Group,
        chr, bin.start
    ) %>% 
    mutate(chr2=chr) %>% 
    # remove entries with < 2 boundaries
    nest(boundaries=c(bin.start)) %>% 
    rowwise() %>% 
    filter(nrow(boundaries) > 1) %>% 
    # convert boundaries to start/end format
    mutate(TADs=list(convert_boundaries_to_TADs(boundaries=boundaries))) %>% 
    ungroup() %>% 
    select(-c(boundaries)) %>% 
    unnest(TADs) %>% 
    mutate(length=end - start) %>% 
    nest(TADs=c(chr, start, end, length)) %>% 
    dplyr::rename(
        'TAD.method'=method,
        'chr'=chr2
    )
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
    # cooltools.TADs.df <- 
    #     load_cooltools_results_for_TADCompare() %>% 
    #     add_column(TAD.method='cooltools')
    # ConsensusTAD TAD results 
    consensusTAD.TADs.df <- 
        load_ConsensusTAD_results_for_TADCompare()
    # default TADCompare method estimates TADs itself, include nothing
    spectralTAD.TADs.df <- 
        expand_grid(
            SampleID=unique(hiTAD.TADs.df$SampleID),
            chr=CHROMOSOMES,
            resolution=unique(hiTAD.TADs.df$resolution)
        ) %>% 
        add_column(
            TADs=NULL, # will be estimated by TADCompare
            TAD.params=NULL,
            TAD.method='spectralTAD'
        )
    # Bind everything together
    bind_rows(
        hiTAD.TADs.df,
        # cooltools.TADs.df,
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
    dplyr::select(-c(TAD.set.index)) %>% 
    dplyr::rename('pre_tads'=TADs)
}

###################################################
# Generate TADCompare results
###################################################
run_TADCompare <- function(
    filepath.Numerator,
    SampleID.Numerator,
    pre_tads.Numerator,
    filepath.Denominator,
    SampleID.Denominator,
    pre_tads.Denominator,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    ...){
    # paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', collapse='; ')
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
    pre_tads <- 
        if (is.null(pre_tads.Numerator)) {
            NULL
        } else {
            list(pre_tads.Numerator, pre_tads.Denominator)
        }
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
    # tad.compare.results$Boundary_Scores %>% as_tibble()
    # tad.compare.results$TAD_Frame %>% as_tibble()
    tad.compare.results$Boundary_Scores %>% 
    as_tibble() %>%
    full_join(
        tad.compare.results$TAD_Frame %>%
        as_tibble() %>% 
        add_column(isTADBoundary=TRUE),
        suffix=c('.All', '.TADs'),
        by=join_by(Boundary)
    ) %>% 
        # {.} -> tcr; tcr
        # tcr %>% count(isTADBoundary, Differential.All, Differential.TADs, Type.All, Type.TADs)
        # tcr %>% 
    mutate(
        isTADBoundary=ifelse(is.na(isTADBoundary), FALSE, isTADBoundary),
        Differential=
            case_when(
                is.na(Differential.TADs) ~ Differential.All,
                TRUE                     ~ Differential.TADs
            ),
        is.Differential=!grepl('Non-Differential', Differential),
        Type=
            case_when(
                is.na(Type.TADs) ~ Type.All,
                TRUE             ~ Type.TADs
            ),
        Enriched.Condition=
            case_when(
                Enriched_In.TADs == 'Matrix 1' ~ SampleID.Numerator,
                Enriched_In.TADs == 'Matrix 2' ~ SampleID.Denominator,
                Enriched_In.All  == 'Matrix 1' ~ SampleID.Numerator,
                Enriched_In.All  == 'Matrix 2' ~ SampleID.Denominator,
                TRUE                      ~ NA
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
    dplyr::rename(
        'TAD.Score.Numerator'=TAD_Score1,
        'TAD.Score.Denominator'=TAD_Score2,
        'TAD.isDifferential'=Differential,
        'TAD.Difference.Type'=Type
    ) %>% 
    rename_with(~ str_replace_all(.x, '_', '.')) %>% 
    dplyr::select(-c(ends_with('.All'), ends_with('.TADs')))
}

run_all_TADCompare <- function(
    comparisons.df,
    hyper.params.df,
    force_redo=FALSE,
    ...){
    # force_redo=TRUE;
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TAD_DIR,
                'results_TADCompare',
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('TAD.method_{TAD.method}'),
                glue('TAD.params_{TAD.params}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                # glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.Numerator}_vs_{SampleID.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    arrange(desc(resolution), desc(chr)) %>% 
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
                TAD.method, 
                TAD.params,
                resolution,
                SampleID.Numerator,
                SampleID.Denominator
            )
        )
    } %>%
        # {.} -> tmp
    # future_pmap(
    pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_TADCompare,
                    # all args also passed as input arguments to run_all*() by pmap
                    ...  # passed from the call to this wrapper()
                )
            },
        ...,  # passed from the call to from run_all_TADCompare
        .progress=TRUE
    )
}

load_TADCompare_results <- function(filepath, ...){
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) # %>%
        # dplyr::count(isTADBoundary, TAD.isDifferential)
    # mutate(
    #     isTADBoundary= 
    #         case_when(
    #             isTADBoundary == 'Differential' ~ TRUE,
    #             isTADBoundary == 'Non-Differential' ~ FALSE,
    #             TRUE ~ isTADBoundary,
    #             # TRUE ~ NA
    #             # isTADBoundary == 'Differential' ~ 'TRUE',
    #             # isTADBoundary == 'Non-Differential' ~ 'FALSE',
    #             # TRUE ~ isTADBoundary
    #         ),
    #     TAD.isDifferential=
    #         case_when(
    #             TAD.isDifferential == 'Differential' ~ TRUE,
    #             TAD.isDifferential == 'Non-Differential' ~ FALSE,
    #             # TRUE ~ TAD.isDifferential
    #         )
    # ) %>%
    #     # dplyr::count(isTADBoundary, TAD.isDifferential)
    # # Only keep TAD boundaries (differential or not) or differential non-boundary bins
    # filter(!is.na(isTADBoundary) & !is.na(TAD.isDifferential)) %>% 
    # filter(!(!isTADBoundary & !TAD.isDifferential))
}

list_all_TADCompare_results <- function(){
    # Get a list of all results files
    file.path(TAD_DIR, 'results_TADCompare') %>%
    parse_results_filelist(
        suffix='-TADCompare.tsv',
        filename.column.name='pair.name'
    ) %>% 
    # Split title into pair of groups ordered by numerator/denominator
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=c('Sample.Group.Numerator', 'Sample.Group.Denominator')
    )
}

load_all_TADCompare_results <- function(...){
    list_all_TADCompare_results() %>% 
    # filter(TAD.method != 'cooltools') %>% 
    # Filter relevant results 
        # {.} -> tmp
        # row_index=2884; tmp$filepath[[row_index]]; filepath <- tmp$filepath[[row_index]]
    mutate(
        results=
            pmap(
            # future_pmap(
                .,
                load_TADCompare_results,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    select(-c(filepath))
}

post_process_TADCompare_results <- function(results.df){
    results.df %>%
    filter(TAD.method != 'cooltools') %>% 
    mutate(
        across(
            c(
                Sample.Group.Numerator,
                Sample.Group.Denominator,
                Enriched.Condition
            ),
            ~ str_remove(.x, '.Merged.Merged')
        ),
        comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}'),
    ) %>% 
    rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=TAD.isDifferential, 
        'Difference'=TAD.Difference.Type
    )
}

