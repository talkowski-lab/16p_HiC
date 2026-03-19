# library(tidyverse)
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
    # boundaries=tmp$boundaries[[1]]; start.col.name='start'; end.col.name='end'
    # convert list of TAD boundaries to start/end format
    # every bnoundary is a start+end excpet for the first and last ones
    tibble(
        TAD.start=boundaries %>% deframe() %>% {.[1:length(.)-1]},
        TAD.end=boundaries %>% deframe() %>% {.[2:length(.)]}
    ) %>%
    add_row(TAD.start=NA, TAD.end=NA) %>% 
    dplyr::rename(
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
    file.path(TAD_RESULTS_DIR, 'method_cooltools')  %>% 
    parse_results_filelist(suffix='-TAD.tsv') %>%
    add_column(method == 'cooltools') %>% 
    # filter(method == 'cooltools') %>% 
    # get_info_from_MatrixIDs(keep_id=FALSE) %>% 
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
    TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='-DI.tsv') %>%
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    filter(method == 'hiTAD') %>% 
    filter(weight == 'balanced') %>%
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

list_all_hiTAD_TADs <- function(){
    TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='-TAD.tsv') %>%
    filter(method == 'hiTAD') %>% 
    get_info_from_MatrixIDs(keep_id=FALSE)
}

load_all_hiTAD_TADs <- function(){
    list_all_hiTAD_TADs() %>% 
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
    filter(weight == 'balanced') %>% 
    filter(isMerged) %>% 
    mutate(length=end - start) %>% 
    # mutate(Sample.Group=str_replace_all(SampleID, '.Merged.Merged', '')) %>% 
    dplyr::select(
        -c(
            weight,
            # threshold,
            # mfvp,
            ReadFilter,
            isMerged,
            Edit, Celltype, Genotype, CloneID, TechRepID,
            # Sample.Group
            # SampleID
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
    # isMerged=tmp$isMerged[[1]]; samples.df=tmp$samples.df[[1]]; Sample.Group=tmp$Sample.Group[[1]]; resolution=tmp$resolution[[1]]; normalization=tmp$normalization[[1]]; z_thresh=tmp$z_thresh[[1]]; window_size=tmp$window_size[[1]]; gap_thresh=tmp$gap_thresh[[1]]; chr=tmp$chr[[1]]; range1=tmp$range1[[1]]; range2=tmp$range2[[1]]; output_dir=tmp$output_dir[[1]]; results_file=tmp$results_file[[1]]
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
    dplyr::rename(
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
                TAD_RESULTS_DIR,
                'method_ConsensusTAD',
                # glue('merged_{isMerged}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                # glue('resolution.type_{resolution.type}'),
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
    arrange(desc(resolution)) %>% 
        # {.} -> tmp; tmp
    # pmap(
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_ConsensusTAD
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_ConsensusTAD,
                    # all extra args also passed as input to run_ConsensusTADs() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this wrapper function by pmap()
        .progress=TRUE
    )
}

load_ConsensusTADs <- function(filepath, ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', sep='', collapse='; ')
    # row_index=586; filepath=tmp$filepath[[row_index]]; method=tmp$method[[row_index]]; z.thresh=tmp$z.thresh[[row_index]]; window.size=tmp$window.size[[row_index]]; gap.thresh=tmp$gap.thresh[[row_index]]; resolution=tmp$resolution[[row_index]]; region=tmp$region[[row_index]]; Sample.Group=tmp$Sample.Group[[row_index]]
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
    ) %>%
    # to uniquely identify rows for pivoting
    mutate(idx=row_number()) %>% 
    mutate(across(ends_with('.score'), as.numeric)) %>% 
    pivot_longer(
        ends_with('.score'),
        names_to='SampleID',
        values_to='score'
    ) %>%
    filter(SampleID == 'consensus.score') %>%
    filter(!is.na(score)) %>% 
    pivot_wider(
        names_from=SampleID,
        values_from=score
    ) %>%
    select(-c(idx))
}

load_all_ConsensusTAD_TADs <- function(...){
    TAD_RESULTS_DIR %>% 
    parse_results_filelist(
        suffix='-ConsensusTADs.tsv',
        filename.column.name='Sample.Group',
    ) %>% 
    filter(method == 'ConsensusTAD') %>% 
        # {.} -> tmp; tmp
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
    # Only keep boundaries
    # mutate(length=end - start) %>% 
    filter(isConsensusBoundary) %>% 
    dplyr::rename('chr'=region)
}

###################################################
# Standardize results across methods
###################################################
load_hiTAD_results_for_TADCompare <- function(force.redo=FALSE){
    check_cached_results(
        results_file=HITAD_TAD_RESULTS_FILE,
        force_redo=force.redo,
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

load_cooltools_results_for_TADCompare <- function(force.redo=FALSE){
    # Load boundary annotations
    check_cached_results(
        results_file=COOLTOOLS_TAD_RESULTS_FILE,
        force_redo=force.redo,
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

load_ConsensusTAD_results_for_TADCompare <- function(force.redo=FALSE){
    check_cached_results(
        results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
        force_redo=force.redo,
        results_fnc=load_all_ConsensusTAD_TADs
    ) %>% 
        # {.} -> tmp; tmp
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

load_all_TAD_results_for_TADCompare <- function(force.redo=FALSE){
    # hiTAD TAD results
    hiTAD.TADs.df <- 
        load_hiTAD_results_for_TADCompare(force.redo=force.redo) %>% 
        add_column(
            TAD.params=NULL,
            TAD.method='hiTAD'
        )
    # cooltools boundary results
    # cooltools.TADs.df <- 
    #     load_cooltools_results_for_TADCompare(force.redo=force.redo) %>% 
    #     add_column(TAD.method='cooltools')
    # ConsensusTAD TAD results 
    consensusTAD.TADs.df <- 
        load_ConsensusTAD_results_for_TADCompare(force.redo=force.redo)
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

load_all_TAD_results <- function(force.redo=FALSE){
    load_all_TAD_results_for_TADCompare(force.redo=force.redo) %>% 
    filter(TAD.method != 'spectralTAD') %>% 
    dplyr::rename('boundaries'=pre_tads) %>%
    select(-c(chr)) %>% 
    unnest(boundaries)
}

