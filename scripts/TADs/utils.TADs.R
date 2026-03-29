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

compute_TAD_stats <- function(
    grouped.df,
    resolution){
    # compute summary stats for each TAD for all bins inside each TAD
    # the score is whatever score is used to decide if a bin is a TAD or not e.g. Insulation, Directionality
    grouped.df %>% 
    # must have columns start, end, bin.start, bin.end, score
    summarize(
        TAD.start.score=sum(ifelse(start == bin.start, bin.score, 0)),
        TAD.end.score=sum(ifelse(end == bin.end, bin.score, 0)),
        across(
            .cols=c(bin.score),
            .fn=
                list(
                    'min'=min,
                    'q25'=partial(stats::quantile, probs=0.25, na.rm=TRUE),
                    'mean'=mean,
                    'median'=median,
                    'q75'=partial(stats::quantile, probs=0.75, na.rm=TRUE),
                    'max'=max,
                    'total'=sum
                ),
            .names="TAD.inner.{.fn}"
        )
    ) %>%
    mutate(
        TAD.length=end - start,
        TAD.bins=TAD.length / resolution
    ) %>% 
    ungroup()
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
load_hiTAD_TADs <- function(
    TAD.filepath,
    DI.filepath,
    resolution,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; '); row.index=1
    # method=tmp$method[[row.index]]; resolution=tmp$resolution[[row.index]]; Edit=tmp$Edit[[row.index]]; Celltype=tmp$Celltype[[row.index]]; Genotype=tmp$Genotype[[row.index]]; CloneID=tmp$CloneID[[row.index]]; TechRepID=tmp$TechRepID[[row.index]]; weight=tmp$weight[[row.index]]; DI.filepath=tmp$DI.filepath[[row.index]]; TAD.filepath=tmp$TAD.filepath[[row.index]]; SampleID=tmp$SampleID[[row.index]]
    # Load bin-wise Adaptibe Directinality scores
    DIs <- 
        read_tsv(
            DI.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'bin.start',
                    'bin.end',
                    'bin.score'
                )
        )
    # load actual TAboundaries
    TADs <- 
        read_tsv(
            TAD.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'start',
                    'end'
                )
        )
    # Annotatte ADI summary stats for each TAD
    TADs %>% 
    # Map scores to TADs
    left_join(
        DIs,
        suffix=c('.TAD', '.DI'),
        by=
            join_by(
                chr,
                within(y$bin.start, y$bin.end, x$start, x$end)
            )
    ) %>% 
    # for each TAD compute summary stats over the bin-wise scores
    group_by(chr, start, end) %>% 
    compute_TAD_stats(resolution=resolution)
}

list_all_hiTAD_TADs <- function(){
    HITAD_TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='.tsv') %>%
    filter(weight == 'balanced') %>% 
    add_column(method='hiTAD') %>% 
    separate_wider_delim(
        MatrixID,
        delim='-',
        names=c('MatrixID', 'feature.type')
    ) %>% 
    pivot_wider(
        names_from=feature.type,
        names_glue="{feature.type}.filepath",
        values_from=filepath
    ) %>% 
    get_info_from_MatrixIDs(
        MatrixID.fields=c('Edit', 'Celltype', 'Genotype', 'CloneID', 'TechRepID', NA, NA, NA)
    )
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
    select(-c(ends_with('.filepath')))
}

post_process_hiTAD_TAD_results <- function(results.df){
    results.df %>% 
    # Subset to only relevant parameters
    filter(weight == 'balanced') %>% 
    dplyr::select(
        -c(
            weight,
            # threshold, mfvp,
            # z.thresh, window.size, gap.thresh,
            Edit, Celltype, Genotype, 
            CloneID, TechRepID,
        )
    )
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

load_ConsensusTAD_TADs <- function(
    filepath,
    resolution,
    ...){
    # row_index=1; paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', sep='', collapse='; ')
    # filepath=tmp$filepath[[row_index]]; z.thresh=tmp$z.thresh[[row_index]]; window.size=tmp$window.size[[row_index]]; gap.thresh=tmp$gap.thresh[[row_index]]; resolution=tmp$resolution[[row_index]]; region=tmp$region[[row_index]]; Edit=tmp$Edit[[row_index]]; Celltype=tmp$Celltype[[row_index]]; Genotype=tmp$Genotype[[row_index]]; Sample.Group=tmp$Sample.Group[[row_index]]; method=tmp$method[[row_index]]
    # Load all the bin-wise + sample-wise scores TAD scores
    scores.df <- 
        filepath %>% 
        read_tsv(
            show_col_types=FALSE,
            progress=FALSE,
        ) %>%
        # to uniquely identify rows (boundaries) for pivoting
        mutate(
            idx=row_number(),
            across(ends_with('.score'), as.numeric)
        ) %>% 
        # has 1 column for every sample's individual TAD score + consensus score
        pivot_longer(
            ends_with('.score'),
            names_to='SampleID',
            values_to='score'
        ) %>%
        # only keep consensus score for each bin
        filter(SampleID == 'consensus.score') %>%
        # filter(!is.na(score)) %>% 
        pivot_wider(
            names_from=SampleID,
            values_from=score
        ) %>%
        dplyr::rename('bin.score'=consensus.score) %>% 
        mutate(bin.end=bin.start + resolution) %>% 
        select(-c(idx))
    # Now get list of boundaries and transform to TAD start:end pairs
    TADs.df <- 
        scores.df %>% 
        mutate(isConsensusBoundary=as.logical(isConsensusBoundary)) %>% 
        filter(isConsensusBoundary) %>% 
        # this is a list of boundaries, transform so each row is a "TAD" i.e. 
        # all contiguous bins between a boundary and the next boundary are the TAD
        mutate(convert_boundaries_to_TADs(boundaries=select(., bin.start))) %>%
        filter(!(is.na(start) & is.na(end))) %>% 
        select(start, end)
    # scores.df %>% nrow()
    # TADs.df %>% nrow()
    # Now map score data for each TAD together and compute summary stats
    TADs.df %>% 
    # Map scores to TADs
    left_join(
        scores.df,
        suffix=c('.TAD', '.bin'),
        by=
            join_by(
                within(y$bin.start, y$bin.end, x$start, x$end)
            )
    ) %>% 
    # for each TAD compute summary stats over the bin-wise scores
    group_by(start, end) %>% 
    compute_TAD_stats(resolution=resolution)
}

list_all_ConsensusTAD_TADs <- function(){
    CONSENSUSTAD_TAD_RESULTS_DIR %>% 
    parse_results_filelist(
        filename.column='Sample.Group',
        suffix='-ConsensusTADs.tsv'
    ) %>%
    add_column(method='ConsensusTAD') # %>% 
    # get_info_from_SampleIDs(
    #     SampleID.col='Sample.Group',
    #     SampleID.fields=
    #         c(
    #             'Edit',
    #             'Celltype',
    #             'Genotype'
    #         )
    # )
}

load_all_ConsensusTAD_TADs <- function(){
    list_all_ConsensusTAD_TADs() %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_ConsensusTAD_TADs,
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
        gap.thresh
    ) %>% 
    # dplyr::select(-c(z.thresh, window.size, gap.thresh)) %>% 
    # Only keep boundaries
    dplyr::rename('chr'=region)
}

###################################################
# Standardize results across methods
###################################################
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

load_all_TAD_results <- function(force_redo_sub=FALSE){
    # Load hiTAD results
    hiTAD.TADs.df <- 
        check_cached_results(
            results_file=HITAD_TAD_RESULTS_FILE,
            force_redo=force_redo_sub,
            # force_redo=TRUE,
            results_fnc=load_all_hiTAD_TADs
        ) %>% 
        post_process_hiTAD_TAD_results() %>%
        dplyr::rename('Sample.Group'=SampleID) %>% 
        add_column(TAD.params=NULL)
    # Load ConsensusTAD results
    consensusTAD.TADs.df <- 
        check_cached_results(
            results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
            force_redo=force_redo_sub,
            # force_redo=TRUE,
            results_fnc=load_all_ConsensusTAD_TADs
        ) %>% 
        post_process_ConsensusTAD_TAD_results()
    # Bind evertything together
    bind_rows(
        hiTAD.TADs.df,
        consensusTAD.TADs.df
    )
}

load_all_TAD_results_for_TADCompare <- function(
    force.redo=FALSE,
    force_redo_sub=FALSE){
    # hiTAD TAD results
    all.TADs.df <- 
        check_cached_results(
            results_file=ALL_TAD_RESULTS_FILE,
            force_redo=force.redo, force_redo_sub=force_redo_sub,
            results_fnc=load_all_TAD_results
        ) %>% 
        select(
            resolution,
            Sample.Group,
            method, TAD.params,
            chr, start, end, TAD.length
        ) %>% 
        mutate(chr.copy=chr) %>% 
        # mutate(length=end - start) %>% 
        nest(TADs=c(chr, start, end, TAD.length)) %>% 
        dplyr::rename(
            'TAD.method'=method,
            'chr'=chr.copy
        )
    # default TADCompare method estimates TADs itself, include nothing
    spectralTAD.TADs.df <- 
        expand_grid(
            SampleID=unique(all.TADs.df$SampleID),
            chr=CHROMOSOMES,
            resolution=unique(all.TADS.df$resolution)
        ) %>% 
        add_column(
            TADs=NULL, # will be estimated by TADCompare
            TAD.params=NULL,
            TAD.method='spectralTAD'
        )
    # Bind everything together
    bind_rows(
        all.TADS.df,
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

pivot_TADs_to_boundaries <- function(results.df){
    # Pivot so the bins that start and end at the boundary position are both included
    results.df %>% 
    dplyr::select(-starts_with('TAD.inner.')) %>% 
    mutate(TAD.index=row_number()) %>% 
    dplyr::rename(
        "start.boundary"=start,
        "end.boundary"=end,
        "start.score"=TAD.start.score,
        "end.score"=TAD.end.score
    ) %>% 
    pivot_longer(
        c(start.boundary, end.boundary, start.score, end.score),
        names_to='stat',
        values_to='value'
    ) %>%
    separate_wider_delim(
        stat,
        delim='.',
        names=c('boundary.side', 'value.type')
    ) %>% 
    pivot_wider(
        names_from=value.type,
        values_from=value
    )
}

