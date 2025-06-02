library(tidyverse)
library(magrittr)
library(tictoc)
library(glue)
library(HiCExperiment)
library(hictkR)
###############
# Utilities
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        message("No results file, will just return data")
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        # Set read/write functions based on filetype
        output_filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (output_filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output_filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {extension}'))
        }
        if (file.exists(results_file) & !force_redo) {
            message(glue('{results_file} exists, not recomputing results'))
            if (return_data) {
                message('Loading cached results...')
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                message(glue('{results_file} exists, recomputing results anyways'))
            } else {
                message(glue('No cached results, generating: {results_file}'))
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    # Or dont save the results, just return the results
    toc()
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

count_dirs <- function(filepath){
    filepath %>%
    dirname() %>% 
    str_split('/') %>% 
    length()
}

parse_results_filelist <- function(
    input_dir,
    suffix,
    filename.column.name,
    param_delim='_',
    ...){
    # !!NOTICE!!
    # This will break if any parameter_dir name has a param_delim character in the name or value
    # This shouldnt  break if the filename has a single param_delim character in it 
    # input_dir=HICREP_DIR; suffix='-hicrep.txt'; filename.column.name='file.pair'; param_delim='_'
    suffix_pattern <- glue('*{suffix}$')
    # List all results files that exist
    input_dir %>% 
    list.files(
        pattern=suffix_pattern,
        recursive=TRUE,
        full.names=FALSE
    ) %>% 
    tibble(fileinfo=.) %>%
    mutate(filepath=file.path(input_dir, fileinfo)) %>% 
    # Extract param info from directory names into structured columns
    separate_longer_delim(
        fileinfo,
        delim='/'
    ) %>%
    mutate(
        fileinfo=
            ifelse(
                grepl(suffix_pattern, fileinfo),
                paste(
                    filename.column.name,
                    fileinfo,
                    sep=param_delim
                ),
                fileinfo
            )
    ) %>% 
    separate_wider_delim(
        fileinfo,
        delim='_',
        too_many="merge",  # in case filenames have delim inside
        names=
            c(
                'Parameter',
                'Value'
            )
    ) %>%
    pivot_wider(
        names_from=Parameter,
        values_from=Value
    ) %>% 
    # Fix column types based on content
    readr::type_convert()
}

process_matrix_name <- function(
    df,
    filename.column.name='matrix.name',
    ...){
    df %>% 
    separate_wider_delim(
        all_of(filename.column.name),
        delim='.',
        too_many='merge',
        names=c(
            'Edit',
            'Genotype',
            'SampleNumber',
            'Celltype',
            NA,
            NA,
            'ReadFilter',
            NA
        )
    ) %>% 
    # Matches Sample.ID in SAMPLE_METADATA_FILE
    mutate(Sample.ID=glue('{Edit}.{Genotype}.{SampleNumber}.{Celltype}')) %>% 
    select(-c(Edit, Genotype, SampleNumber, Celltype))
}
###############
# Annotate contacts with specified regions
annotate_contact_region_pairs <- function(
    contacts.df,
    regions.of.interest,
    most_specific_only=TRUE,
    ...){
    # regions of interest is a tibble where each row is a specific contiguous genomic range
    # see GENOMIC_REGIONS variable for example
    # For each bin (1, left and 2, right) in each pair of bins annotate the region is belongs to
    region.pairs.of.interest <- 
        regions.of.interest %>%
        add_column(join_dummy=0) %>% 
        # all possible pairs of regions to annotate
        # each bin-pair will be assigned 1 region pair, even if regions overlap/contain eachother
        full_join(
            .,
            {.},
            suffix=c('1','2'),
            relationship='many-to-many',
            by='join_dummy'
        ) %>% 
        filter(region.chr1 == region.chr2) %>% 
        select(-c(join_dummy)) %>% 
        # factor order so most specific region interactions have the lowest factor level
        mutate(
            region=
                pmap_chr(
                    list(region1, region2),
                    ~ paste(sort(c(...)), collapse=' vs ')
                ),
            region.dist=region.dist1 + region.dist2,
            region=fct_reorder(region, region.dist) 
        )
    # contacts.df should be a tibble where each row is a pair of genomic bins i.e. output of
    # the load_mcool_file(s) functions
    contacts.df %>% 
    # Now annotate all region-pairs that each bin-pair belongs to 
    # (bin1 in region1 and bin2 in region2)
    full_join(
        region.pairs.of.interest,
        relationship='many-to-many',
        by=
            join_by(
                chr == region.chr1,
                between(
                    range1,
                    region.start1,
                    region.end1
                ),
                between(
                    range2,
                    region.start2,
                    region.end2
                )
            )
    ) %>% 
    # Pick most specific (smallest) annotated region for each bin if specified
    {
        if (most_specific_only) {
            group_by(
                .,
                Sample.ID,
                chr,
                range1,
                range2
            ) %>% 
            slice_min(
                region.dist,
                n=1,
                with_ties=FALSE
            ) %>% 
            ungroup()
        } else {
            .
        }
    } %>% 
    # Clean up 
    select(-c(starts_with('region.')))
}

load_annotated_contacts_pairs <- function(
    pattern,
    range1s,
    resolutions,
    regions.of.interest,
    ...){
    # Load all contacts
    load_mcool_files(
        pattern=pattern,
        range1s=range1s,
        resolutions=resolutions,
    ) %>%
    mutate(is.Merged=grepl('Merged', Sample.ID)) %>%
    separate_wider_delim(
        Sample.ID,
        delim=fixed('.'),
        cols_remove=FALSE,
        names=c(
            'Edit',
            'Genotype',
            'SampleNumber',
            'Celltype'
        )
    ) %>% 
    annotate_contact_regions_pairs(
        regions.of.interest=regions.of.interest,
        most_specific_only=TRUE
    )
}
###############
# Load Data
load_sample_metadata <- function(
    sample_metadata_file=SAMPLE_METADATA_FILE,
    clean=TRUE){
    sample_metadata_file %>%
    read_tsv() %>%
    {
        if (clean) {
            select(
                ., 
                -c(
                    # Sequencing.Run,
                    Lane,
                    Investigator.ID,
                    Sample.Info,
                    Fastq.R1,
                    Fastq.R2
                    # Date,
                    # Edit,
                    # Genotype,
                    # Sample.Number,
                    # Celltype,
                    # Sample.ID,
                    # Sample.Name,
                )
            )

        } else {
            .
        }
    }
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(col_names=c('Chr', 'chr.total.bp'))
}

load_rgd_regions <- function(
    normalizations=c("NONE", "ICE"),
    window.size=1e6,
    ...){
    GENOMIC_REGIONS %>% 
    filter(region.group == "RGDs") %>% 
    # context window of 1Mb around each region
    add_column(join_dummy=0) %>% 
    full_join(
        .,
        expand_grid(
            join_dummy=0,
            window.size=window.size,
            normalization=normalizations,
            cis=TRUE
        ),
        relationship='many-to-many'
    ) %>% 
    # format querys for loading contacts via fetch()
    rowwise() %>% 
    mutate(
        range1=glue('{region.chr}:{max(0, region.start - window.size)}-{region.end + window.size}'),
        range2=range1
    ) %>%
    ungroup() %>% 
    select(-c(region.group, join_dummy))
}

load_rgd_contacts <- function(
    region.df,
    resolutions,
    ...){
    load_mcool_files(
        pattern='*.mapq_30.1000.mcool',
        resolutions=resolutions,
        region.df=region.df,
        range1s=NULL,
        range2s=NULL,
        progress=TRUE
    ) %>% 
    # group_by(Sample.ID, region) %>% 
    # add_tally(wt=IF, name='region.total.contacts') %>% 
    # ungroup() %>%
    # group_by(Sample.ID, region) %>% 
    # add_count(name='region.nbins') %>% 
    # ungroup() %>%
    select(-c(region.chr, weight, cis))
}
###############
# Load Files
load_genome_coverage <- function(
    weights=NULL,
    resolutions=NULL,
    ...){
    # List input files (generated by cooltools coverage)
    parse_results_filelist(
        input_dir=COVERAGE_DIR,
        suffix='-coverage.tsv',
        filename.column.name='matrix.name',
        param_delim='_',
    ) %>%
    process_matrix_name() %>% 
    # Ignore results not matching param values
    { 
        if (is.null(resolutions)) { 
            . 
        } else { 
            filter(., resolution %in% resolutions) 
        } 
    } %>%
    { 
        if (is.null(weights)) { 
            . 
        } else { 
            filter(., weight %in% weights) 
        } 
    } %>%
    mutate(
        coverage=
            purrr::pmap(
                .l=.,
                .f=function(filepath, ...){
                    read_tsv(
                        filepath,
                        show_col_types=FALSE,
                        progress=FALSE
                    )
                }
            )
    ) %>%
    unnest(coverage) %>%
    rename_with(~ str_remove(.x, '_(raw|weight)')) %>% 
    rename(
        'chr'=chrom,
        'coverage.cis'=cov_cis,
        'coverage.total'=cov_tot,
    ) %>% 
    mutate(chr=factor(chr, levels=CHROMOSOMES)) %>%  
    pivot_longer(
        starts_with('coverage'),
        names_to='metric',
        names_prefix='coverage.',
        values_to='coverage'
    )
}

load_mcool_file <- function(
    filepath,
    resolution,
    normalization,
    range1="",
    range2="",
    cis=TRUE,
    ...){
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        normalization=normalization,
        join=TRUE,
        query_type='UCSC',
        type='df'
    ) %>% 
    as_tibble() %>%
    add_column(weight=normalization) %>% 
    # format column names
    {
        if (cis) {
            filter(., chrom1 == chrom2) %>% 
            select(
                -c(
                    chrom2,
                    end1,
                    end2
                )
            ) %>% 
            rename(
                'chr'=chrom1,
                'range1'=start1,
                'range2'=start2,
                'IF'=count
            )
        } else {
            rename(
                'chr1'=chrom1,
                'chr2'=chrom2,
                'range1'=start1,
                'range2'=start2,
                'IF'=count
            )
        }
    }
}

load_mcool_files <- function(
    pattern,
    resolutions,
    region.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    ...){
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    mutate(matrix.name=basename(filepath)) %>% 
    mutate(join_dummy=0) %>% 
    {
        if (is.null(region.df)) {
            if (is.null(range2s)) {
                full_join(
                    expand_grid(
                        .,
                        resolution=resolutions,
                        range1=range1s,
                        join_dummy=0
                    ) %>%
                    mutate(range2=range1),
                    relationship='many-to-many',
                )
            } else {
                full_join(
                    .,
                    expand_grid(
                        resolution=resolutions,
                        range1=range1s,
                        range2=range2s,
                        join_dummy=0
                    ),
                    relationship='many-to-many',
                )
            }
        } else {
            full_join(
                .,
                region.df %>% 
                add_column(join_dummy=0) %>%
                full_join(
                    expand_grid(
                        resolution=resolutions,
                        join_dummy=0
                    ), 
                    by=join_by(join_dummy)
                )
            )
        }
    } %>% 
    process_matrix_name() %>% 
    mutate(resolution=as.integer(resolution)) %>% 
    mutate(
        contacts=
            purrr::pmap(
                .l=.,
                .f=load_mcool_file,
                .progress=progress
            )
    ) %>% 
    select(
        -c(
            filepath,
            join_dummy,
            range1,
            range2
        )
    ) %>%
    unnest(contacts)
}

load_sparse_matrix_file <- function(
    filename,
    MHC_format=FALSE){
    # input format compatible with multiHiCCompare functions
    if (MHC_format) {
        filename %>% 
        read.table(
            .,
            header=FALSE
        ) %>% 
        {.[, c(1, 2, 5, 7)]} %>% 
        setNames(
            c(
                'chr', 
                'range1',
                'range2',
                'IF'
            )
        )
    } else {
        filename %>% 
        read.table(
            header=FALSE,
            col.names=
                c(
                  'A.chr',
                  'A.start',
                  'A.end',
                  'B.chr',
                  'B.start',
                  'B.end',
                  'IF'
                )
        )
    }
}

load_sparse_matrix_files <- function(
    input_dir,
    suffix='-sparse.matrix.txt',
    param_names,
    MHC_format,
    ...){
    input_dir %>% 
    parse_results_filelist(
    )
}
