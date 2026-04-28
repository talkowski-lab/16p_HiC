###################################################
# Dependencies
###################################################
# library(tidyverse)
library(stringi)
library(glue)
library(furrr)
library(TADCompare)

###################################################
# Utilities
###################################################
convert_boundaries_to_TADs <- function(
    bins.df,
    include.tails=FALSE,
    start.col.name='start',
    end.col.name='end',
    ...){
    # bins.df=tmp2$bins.df[[1]]
    # convert list of TAD boundaries to start/end format
    # every boundary is a start+end excpet for the first and last ones
    bins.df %>% 
    filter(is.boundary) %>%
    pull(bin.start) %>%
    {
        if (include.tails) {
            c(
                bins.df %>% pull(bin.start) %>% min(), # first bin on chr
                .,                                     # all annotated boundaries
                bins.df %>% pull(bin.start) %>% max()  # last bin on chr
            )
        } else {
            .
        }
    } %>% 
    tibble(boundaries=.) %>%
    rename(!!start.col.name := boundaries) %>% 
    mutate(!!end.col.name := lead(!!sym(start.col.name))) %>% 
    filter(!is.na(start), !is.na(end))
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
            .fns=
                list(
                    'min'=   partial(min,    na.rm=TRUE),
                    'q25'=   partial(stats::quantile, probs=0.25, na.rm=TRUE),
                    'mean'=  partial(mean,   na.rm=TRUE),
                    'median'=partial(median, na.rm=TRUE),
                    'q75'=   partial(stats::quantile, probs=0.75, na.rm=TRUE),
                    'max'=   partial(max,    na.rm=TRUE),
                    'total'= partial(sum,    na.rm=TRUE)
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

calculate_MoC <- function(
    TADs.P1,
    TADs.P2,
###################################################
# Cooltools 
###################################################
load_cooltools_file <- function(
    filepath,
    boundaries_only=FALSE,
    ...){
    # TADs.P1, TADs.P2 are both tibbles with 
    # the following 3 columns: start, end, length
    # add indices to track all pairs of TADs
    # MoC normalization constant
    nTADs.P1 <- nrow(TADs.P1)
    nTADs.P2 <- nrow(TADs.P2)
    norm_const <- 1 / (sqrt(nTADs.P1 * nTADs.P2) - 1)
    # nTADs.P1; nTADs.P2; norm_const;
    # Now find all overlapping pairs of TADs
    inner_join(
        TADs.P1 %>% mutate(idx=row_number()),
        TADs.P2 %>% mutate(idx=row_number()),
        suffix=c('.P1', '.P2'),
        by=join_by(overlaps(x$start, x$end, y$start, y$end))
    ) %>% 
    # https://link.springer.com/article/10.1186/s13059-018-1596-9#Sec9
    # "Assessment of TAD calller performance"
    rowwise() %>% 
    mutate( 
        # F_ij^2 / (P_i * Q_j), 1 pair of TADs per row
        intersection=min(end.P1, end.P2) - max(start.P1, start.P2),
        moc.inner=((intersection**2) / (TAD.length.P1 * TAD.length.P2))
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
    # group_by(idx.P1) %>% slice_max(moc.inner) %>%  ungroup() %>% 
    # group_by(idx.P2) %>% slice_max(moc.inner) %>% 
    # When multiple TADs overlap, count only the most overlapping match
    # group_by(idx.P1) %>%
    # slice_max(moc.inner)
    ungroup() %>% 
    summarize(
        n.Overlaps=n(),
        n.TADs.P1=length(unique(idx.P1)),
        n.TADs.P2=length(unique(idx.P2)),
        MoC=(sum(moc.inner) - 1) / (sqrt(n.TADs.P1 * n.TADs.P2) - 1)
        # MoC=(sum(moc.inner) - 1) * norm_const
    )
}

calculate_all_MoCs <- function(
    nested.TADs.df,
    suffixes=NULL,
    delim='.',
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; ')
    # suffixes=c('Numerator', 'Denominator'); delim='.'
    nested.TADs.df %>% 
    enumerate_pairwise_comparisons(
        delim=delim,
        suffixes=c('P1', 'P2'),
        ...
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
    # Finally compute all MoCs for all listed pairs of annotations
    mutate(
        MoCs=
            # pmap(
            future_pmap(
                .l=.,
                .f=calculate_MoC,
                .progress=TRUE
            )
    pivot_wider(
        names_from=stat,
        values_from=value
    ) %>%
    select(-c(TADs.P1, TADs.P2)) %>%
    unnest(MoCs) %>%
    rename('is.boundary'=is_boundary) %>% 
    mutate(is.boundary=as.logical(is.boundary)) %>% 
    {
        if (!is.null(suffixes) & length(suffixes) == 2) {
            rename_with(
                ., 
                .cols=ends_with('P1'),
                ~str_replace(.x, 'P1$', suffixes[[1]])
            ) %>% 
            rename_with(
                .cols=ends_with('P2'),
                ~str_replace(.x, 'P2$', suffixes[[2]])
            )
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
    add_column(TAD.metric='ADI') %>% 
    # for each TAD compute summary stats over the bin-wise scores
    group_by(TAD.metric, chr, start, end) %>% 
    compute_TAD_stats(resolution=resolution)
}

list_all_hiTAD_TADs <- function(){
    HITAD_TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='.tsv') %>%
    # hiTAD only expects balanced matrices as input 
    filter(weight == 'balanced') %>% 
    add_column(TAD.method='hiTAD') %>% 
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
    get_info_from_MatrixIDs(include_fields=FALSE) %>%
    select(-c(SampleID))
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
    dplyr::select(-c(weight))
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
        add_column(TAD.params=NULL)
    # Load cooltools TAD results
    cooltools.TADs.df <- 
        check_cached_results(
            results_file=COOLTOOLS_TAD_RESULTS_FILE,
            force_redo=force_redo_sub,
            # force_redo=TRUE,
            results_fnc=load_all_cooltools_TADs
        ) %>% 
        post_process_cooltools_TAD_results()
    # Load ConsensusTAD results
    # consensusTAD.TADs.df <- 
    #     check_cached_results(
    #         results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
    #         force_redo=force_redo_sub,
    #         # force_redo=TRUE,
    #         results_fnc=load_all_ConsensusTAD_TADs
    #     ) %>% 
    #     post_process_ConsensusTAD_TAD_results()
    # Bind evertything together
    bind_rows(
        hiTAD.TADs.df,
        cooltools.TADs.df
        # consensusTAD.TADs.df
    )
}

post_process_all_TAD_results <- function(results.df){
    results.df %>% 
    filter(TAD.length < 10 * 1e6) %>% 
    mutate(
        TAD.size.band=
            case_when(
                length > 5e6 ~ '>  5Mb',
                length > 2e6 ~ '>  2Mb',
                length > 1e6 ~ '>  1Mb',
                length > 5e5 ~ '>  500Kb',
                TRUE         ~ '<= 500Kb'
            ) %>%
            factor(levels=c('<= 500Kb', '>  500Kb', '>  1Mb', '>  2Mb', '>  5Mb'))
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
            'chr'=chr.copy,
            'TAD.method'=method
        )
    # default TADCompare method estimates TADs itself, include nothing
    spectralTAD.TADs.df <- 
        expand_grid(
            Sample.Group=unique(all.TADs.df$Sample.Group),
            chr=CHROMOSOMES,
            resolution=unique(all.TADs.df$resolution)
        ) %>% 
        add_column(
            TADs=NULL, # will be estimated by TADCompare
            TAD.params=NULL,
            TAD.method='spectralTAD'
        )
    # Bind everything together
    bind_rows(
        all.TADs.df,
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
    ) %>% 
    dplyr::rename(
        "boundary.start"=boundary,
        "boundary.score"=score
    )
}

###################################################
# Generate TADCompare results
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

run_TADCompare <- function(
    filepath.Numerator,
    Sample.Group.Numerator,
    pre_tads.Numerator,
    filepath.Denominator,
    Sample.Group.Denominator,
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
                Enriched_In.TADs == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.TADs == 'Matrix 2' ~ Sample.Group.Denominator,
                Enriched_In.All  == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.All  == 'Matrix 2' ~ Sample.Group.Denominator,
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
                TADCOMPARE_DIR,
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
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    # filter(chr != 'chrY') %>% 
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
                # z_thresh,
                # window_size,
                # gap_thresh,
                # TAD.params,
                TAD.method, 
                resolution,
                Sample.Group.Numerator,
                Sample.Group.Denominator
            )
        )
    } %>%
        # {.} -> tmp; tmp
        # tmp %>% 
    # pmap(
    future_pmap(
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

###################################################
# Load TADCompare results
###################################################
load_TADCompare_results <- function(
    filepath,
    boundaries_only=TRUE,
    ...){
    # row_index=1; filepath=tmp$filepaths[[row_index]][[1]];
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
        # {.} -> ltr.tmp; ltr.tmp
        # ltr.tmp %>% count(isTADBoundary, is.Differential, TAD.Difference.Type)
    {
        if (boundaries_only) {
            filter(., isTADBoundary)
        } else {
            .
        }
    } %>% 
    filter(!is.na(TAD.Difference.Type))
}

load_and_correct_TADCompare_results <- function(
    filepaths,
    nom.threshold,
    # fdr.threshold,
    gw.fdr.threshold,
    ...){
    # filepaths=tmp$filepaths[[1]]
    filepaths %>% 
    mutate(
        results=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                id='tmpID',
                show_col_types=FALSE,
                progress=FALSE
            )
    ) %>% 
    unnest(results) %>% 
    # The score provided by TADCompare is functionallt a z-score distributed at N(0,1)
    # so we can compute a regular p-value to decide if a TAD's boundary score is significantly different
    # between conditions
    # See Section 2.7 here
    # https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00158/full 
    # calculate pvalue from Gap Score (a z-score) calcualted by TADCompare
    mutate(p.value=2 * pnorm(abs(Gap.Score), lower.tail=FALSE)) %>% 
    mutate(p.adj.gw=p.adjust(p.value, method='BH')) %>% 
    ungroup() %>% 
    # only keep sufficiently sifnificant siginicant differences
    filter(
        p.adj.gw < gw.fdr.threshold,
        # p.adj    < fdr.threshold,
        p.value  < nom.threshold
    ) %>%
    select(-c(tmpID))
}

list_all_TADCompare_results <- function(){
    # Get a list of all results files
    TADCOMPARE_DIR %>% 
    parse_results_filelist(
        suffix='-TADCompare.tsv',
        filename.column.name='pair.name'
    ) %>% 
    # Split title into pair of groups ordered by numerator/denominator
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=c('SampleID.Numerator', 'SampleID.Denominator')
    ) %>% 
        # {.} -> tmp; tmp
        # tmp %>% 
    extract_all_sample_pair_metadata(
        SampleID.cols=c('SampleID.Numerator', 'SampleID.Denominator'),
        SampleID.fields=c('Edit', 'Celltype', 'Genotype'),
        suffixes=c('Numerator', 'Denominator')
    )
}

load_all_TADCompare_results <- function(
    nom.threshold,
    # fdr.threshold,
    gw.fdr.threshold,
    ...){
    # gw.fdr.threshold=1; fdr.threshold=0.1; nom.threshold=0.05
    list_all_TADCompare_results() %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
    mutate(
        results=
            # pmap(
            future_pmap(
                .l=.,
                # load_TADCompare_results,
                .f=load_and_correct_TADCompare_results,
                nom.threshold=nom.threshold,
                # fdr.threshold=fdr.threshold,
                gw.fdr.threshold=gw.fdr.threshold,
                boundaries_only=TRUE,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

load_correct_count_TADCompare_results <- function(
    filepaths,
    sig.colname='p.adj.gw',
    ...){
    filepaths %>% 
    load_and_correct_TADCompare_results(
        nom.threshold=1,
        # fdr.threshold=1,
        gw.fdr.threshold=1
    ) %>% 
    # for each thresh, make binary col if TAD difference meets threshold
    mutate(
        "sig.lvl.{sig.colname} < 1e-15" := .data[[sig.colname]] <  1e-15,
        "sig.lvl.{sig.colname} < 1e-10" := .data[[sig.colname]] <  1e-10,
        "sig.lvl.{sig.colname} < 1e-05" := .data[[sig.colname]] <  1e-5,
        "sig.lvl.{sig.colname} < 0.001" := .data[[sig.colname]] <  1e-3,
        "sig.lvl.{sig.colname} < 0.05 " := .data[[sig.colname]] <  0.05,
        "sig.lvl.{sig.colname} < 0.1  " := .data[[sig.colname]] <  0.10,
        "sig.lvl.N.S."                  := .data[[sig.colname]] >= 0.10
        # "sig.lvl.NA"                    := is.na(.data[[sig.colname]])
    ) %>% 
    pivot_longer(
        starts_with('sig.lvl.'),
        names_to='sig.lvl',
        names_prefix='sig.lvl.',
        values_to='meet.sig.lvl'
    ) %>% 
    # Inclusively count how many TAD differences meet each thrshold across categories
    # This produces inclusive counts  for each significance threshold i.e. 
    # the number of TAD differences < 0.1 also includes all differences <= 0.01
    filter(meet.sig.lvl) %>% 
    count(
        isTADBoundary,
        is.Differential,
        TAD.Difference.Type,
        Enriched.Condition,
        region,
        sig.lvl
    )
}

load_correct_count_all_TADCompare_results <- function(){
    list_all_TADCompare_results() %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
    mutate(
        results=
            future_pmap(
                .l=.,
                .f=load_correct_count_TADCompare_results,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

post_process_TADCompare_results <- function(results.df){
    results.df %>%
    # filter(TAD.method != 'cooltools') %>% 
    filter(isBoundary) %>% 
    mutate(
        # isBoundary=ifelse(isBoundary, 'TAD', 'Not TAD'),
        across(
            c(
                SampleID.Numerator,
                SampleID.Denominator,
                Enriched.Condition
            ),
            ~ str_remove(.x, '.Merged.Merged')
        ),
    ) %>% 
    # mutate(log.p.adj.gw=-log10(p.adj.gw)) %>% 
    select(
        -c(
            isDifferential,
            isBoundary,
            z.thresh,
            window.size,
            gap.thresh
        )
    ) %>% 
    relocate(
        c(
            resolution,
            TAD.method,
            TAD.params,
            SampleID.Numerator, SampleID.Denominator, 
            chr,
            DifferenceType,
            Enriched.Condition
        )
    )
}

