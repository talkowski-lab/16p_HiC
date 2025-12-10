library('FitHiC')

###################################################
# Generate FitHiC2 commands to run
###################################################
format_fragments <- function(
    coverage.filepath,
    resolution,
    contact.type,
    ...){
    coverage.filepath %>%
    load_genome_coverage() %>% 
    filter(metric == contact.type) %>% 
    mutate(
        # mappable=,
        chr=rename_chrs(chr, to_numeric=TRUE),
        fragmentMid=bin.start + resolution / 2,
        marginalizedContactCount=coverage
    ) %>%
    add_column(
        extraField=NA,
        mappable=1
    ) %>% 
    select(
        chr,
        extraField,
        fragmentMid,
        marginalizedContactCount,
        mappable
    ) %>% 
    relocate(
        chr,
        extraField,
        fragmentMid,
        marginalizedContactCount,
        mappable
    )
}

format_interactions <- function(
    matrix.filepath,
    resolution,
    normalization,
    chr,
    cis,
    ...){
    matrix.filepath %>% 
    load_mcool_file(
        resolution=resolution,
        normalization=normalization,
        range1=chr,
        cis,
        type='df'
    ) %>% 
    rename(
        'chr1'=chr.,
        'chr2'=chr.,
        'fragmentMid1'=bin.start1,
        'fragmentMid2'=bin.start2,
        'contactCount'=value
    ) %>%
    mutate(
        fragmentMid1=fragmentMid1 + resolution / 2,
        fragmentMid2=fragmentMid2 + resolution / 2
    )
}

run_FitHiC <- function(
    matrix.filepath,
    coverage.filepath,
    bias.filepath,
    Sample.Group,
    chr,
    contact.type,
    output_dir,
    resolution,
    normalization,
    contact_type,
    noOfPasses,
    nOfBins,
    mappabilityThreshold,
    lowerbound,
    upperbound,
    biasLowerBound,
    biasUpperBound,
    # distUpThres,
    # distLowThres,
    ...){ 
    fragments.df <- 
        format_fragments(
            coverage.filepath
            resolution=resolution,
            contact.type=contact.type,
        )
    interactions.df <- 
        format_interactions(
            matrix.filepath,
            resolution=resolution,
            normalization=normalization,
            chr=chr,
            cis=(contact_type == 'cis')
        )
    # cmd <- glue("fithic --passes {noOfPasses} --nOfBins {nOfBins} --mappabilityThres {mappabilityThreshold} --outdir {output_dir} --lib {Sample.Group} --contactType {contact_type} --lowerbound {lowerbound} --upperbound {upperbound} --resolution {resolution} --biasLowerBound {biasLowerBound} --biasUpperBound {biasUpperBound} --biases {bias.filepath} --fragments {fragment.filepath} --interactions {matrix.filepath}")
    FitHiC(
        fragment.filepath,
        matrix.filepath,
        outdir=output_dir,
        biasfile="none",
        libname=Sample.Group,
        noOfPasses=passes,
        nOfBins=nOfBins,
        mappabilityThreshold=mappabilityThreshold,
        distUpThres=distUpThres,
        distLowThres=distLowThres,
        visual=FALSE,
        useHiCPro=FALSE
    )
}

run_all_FitHiC <- function(
    comparisons.df,
    hyper.params.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    ...){
    # chromosomes=CHROMOSOMES; force_redo=TRUE;
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        bias.filepath=,
        fragment.filepath=,
        output_dir=
            file.path(
                LOOP_DIR,
                # glue('merged_{isMerged}'),
                glue('method_{method}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('normalization_{normalization}'),
                glue('passes_{passes}'),
                glue('nOfBins_{nOfBins}'),
                glue('biasLowerBound_{biasLowerBound}'),
                glue('biasUpperBound_{biasUpperBound}'),
                glue('mappabilityThres_{mappabilityThres}'),
                glue('upperbound_{upperbound}'),
                glue('lowerbound_{lowerbound}'),
                glue('contact.type_{contact_type}'),
                glue('region_{chr}'),
                glue('Sample.Group_{Sample.Group}')
            )
    ) %>% 
    arrange(resolution) %>% 
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
                    results_fnc=run_FitHiC,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

load_FitHiC_results <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE,
    )
}

load_all_FitHiC_results <- function(
    ...){
    LOOP_DIR %>%
    parse_results_filelist() %>% 
    get_info_from_MatrixIDs(keep_id=FALSE) %>% 
    mutate(
        Loops=
            # pmap(
            future_pmap(
                .,
                load_FitHiC_results,
                .progress=TRUE
            )
    ) %>%
    unnest(Loops) %>% 
    select(-c(filepath))
}

post_process_FitHiC_results <- function(results.df){
    results.df %>%
    standardize_data_cols()
}

