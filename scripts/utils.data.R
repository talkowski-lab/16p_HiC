###################################################
# Dependencies
###################################################
library(tidyverse)
library(magrittr)
library(tictoc)
library(glue)
library(HiCExperiment)
library(hictkR)

###################################################
# Pairsing and Caching
###################################################
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

parse_results_filelist <- function(
    input_dir,
    suffix,
    filename.column.name='MatrixID',
    pattern=NA,
    param_delim='_',
    ...){
    # !!NOTICE!!
    # This will break if any parameter_dir name has a param_delim character in the name or value, not as the delimiter
    # This shouldnt  break if the filename has a single param_delim character in it 
    # i.e. min_resolution-0.45 is fine when param_delim='-'
    #      min_resolution_0.45 is will break with any value for param_delim
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
    {
        if (is.na(pattern)) {
            .
        } else {
            filter(., grepl(pattern, filepath))
        }
    } %>% 
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

glue_sample_ID <- function(
    Edit,
    Celltype,
    Genotype,
    CloneID,
    TechRepID,
    prefix='',
    ...) {
    glue('{prefix}{Edit}.{prefix}{Celltype}.{prefix}{Genotype}.{prefix}{CloneID}.{prefix}{TechRepID}')
}

get_info_from_SampleIDs <- function(
    df,
    sample_ID_col='SampleID',
    col_prefix='',
    keep_id=TRUE,
    nest_col=NA,
    SampleID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID'
        ),
    ...){
    df %>%
    mutate(
        !!as.character(glue('{col_prefix}isMerged')) :=
            ifelse(
                grepl('Merged', !!sym(sample_ID_col)),
                'Merged',
                'Individual'
            ) %>% 
            factor()
    ) %>% 
    separate_wider_delim(
        all_of(sample_ID_col),
        delim=fixed('.'),
        names=SampleID.fields %>% paste0(col_prefix, .),
        cols_remove=!keep_id
    ) %>% 
    {
        if (!is.na(nest_col)) {
            if (col_prefix  != '') {
                nest(., !!nest_col := starts_with(col_prefix))
            } else {
                nest(
                    ., 
                    !!nest_col :=
                        all_of(
                            c(
                                SampleID.fields,
                                'isMerged'
                            )
                        )
                )
            }
        } else {
            .
        }
    }
}

get_info_from_MatrixIDs <- function(
    df,
    matrix_ID_col='MatrixID',
    sample_ID_col='SampleID',
    col_prefix='',
    keep_id=TRUE,
    nest_col=NA,
    delim=fixed('.'),
    separate_names=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID',
            NA,
            'ReadFilter',
            NA
        ),
    ...){
    # matrix_ID_col='MatrixID.P1'; sample_ID_col='SampleID.P1'; col_prefix='SampleInfo.P1.'; keep_id=FALSE; nest_col=NA
    # Set up some column names/variables
    col_names <- c(separate_names[!is.na(separate_names)], 'isMerged')
    # extract + format metadata
    df %>%
    mutate(isMerged=grepl('Merged', !!sym(matrix_ID_col))) %>% 
    separate_wider_delim(
        all_of(matrix_ID_col),
        delim=delim,
        names=separate_names,
        too_many='merge',
        cols_remove=!keep_id
    ) %>% 
    # create SampleID
    {
        if (!is.null(sample_ID_col)) {
            mutate( ., !!sample_ID_col := pmap(.l=., .f=glue_sample_ID)) %>%
            unnest(!!sym(sample_ID_col))
        } else {
            .
        }
    } %>% 
    # rename columns with a common prefix
    {
        if (col_prefix != '') {
            rename_with(
                ., 
                .fn=~ paste0(col_prefix, .x), 
                .cols=all_of(col_names)
            )
        } else {
            .
        }
    } %>% 
    # nest and rename columns as specified
    {
        if (!is.na(nest_col) && col_prefix != '') {
            nest(
                .,
                !!nest_col := all_of(paste0(col_prefix, col_names))
            )
        } else if (!is.na(nest_col) && col_prefix == '') {
            nest(
                .,
                !!nest_col := all_of(col_names)
            )
        } else {
            .
        }
    }
}

###################################################
# Format stuff
###################################################
scale_numbers <- function(
    numbers,
    accuracy=2,
    force_numeric=FALSE){
    if (is.character(numbers) | is.factor(numbers)) {
        numbers %>%
        as.character() %>% 
        tibble(resolution.str=.) %>%
        mutate(
            suffix=str_extract(resolution.str, '[KMG]b'),
            magnitude=
                case_when(
                    suffix == 'Kb' ~ 1e3,
                    suffix == 'Mb' ~ 1e6,
                    suffix == 'Gb' ~ 1e9,
                    TRUE ~ 1
                ),
            resolution=
                resolution.str %>% 
                str_remove('[KMG]b') %>%
                as.integer() %>%
                multiply_by(magnitude)
        ) %>% 
        pull(resolution)
    } else if (is.numeric(numbers) & !force_numeric) {
        numbers %>%
        tibble(resolution=.) %>%
        mutate(
            magnitude=resolution %>% log10() %>% floor() %>% {. %/% 3} %>% {. * 3},
            suffix=
                case_when(
                    magnitude == 3 ~ 'Kb',
                    magnitude == 6 ~ 'Mb',
                    magnitude == 9 ~ 'Gb',
                    TRUE ~ ''
                ),
            resolution.digits=signif(resolution / 10**magnitude, digits=accuracy),
            resolution.str=glue('{resolution.digits}{suffix}')
        ) %>% 
        mutate(resolution.str=fct_reorder(resolution.str, resolution)) %>% 
        pull(resolution.str)
    } else {
        numbers
    }
}

rename_chrs <- function(chrs){
    case_when(
        chrs == 23           ~ 'X',
        chrs == 24           ~ 'Y',
        chrs > 0 & chrs < 23 ~ as.character(chrs),
        TRUE ~ NA
    ) %>% 
    paste0('chr', .) %>% 
    factor(levels=CHROMOSOMES)
}

standardize_data_cols <- function(results.df){
    # results.df <- HITAD_DI_RESULTS_FILE %>% read_tsv(show_col_types=FALSE)
    # results.df %>% group_by(resolution, weight, SampleID, chr) %>% slice_head(n=1)
    results.df %>% 
    {
        if ('isMerged' %in% colnames(.)) {
            if (is.logical(results.df$isMerged)) {
                mutate(
                    .,
                    isMerged=
                        ifelse(isMerged, 'Merged', 'Individual') %>%
                        factor(levels=c('Merged', 'Individual'))
                )
            } else {
                mutate(., isMerged=factor(isMerged, levels=c('Merged', 'Individual')))
            }
        } else {
            .
        }
    } %>% 
    {
        if ('resolution' %in% colnames(.)) {
            mutate(
                .,
                resolution=
                    resolution %>%
                    scale_numbers(force_numeric=TRUE) %>%
                    scale_numbers(),
            )
        } else {
            .
        }
    } %>% 
    {
        if ('chr' %in% colnames(.)) {
            mutate(., chr=factor(chr, levels=CHROMOSOMES))
        } else {
            .
        }
    }
}

###################################################
# Handle pairs of samples
###################################################
join_all_rows <- function(
    df1,
    df2=NULL,
    join_keys=c(),
    ...){
    if (is_tibble(df2)) {
        full_join(
            df1 %>% add_column(join_dummy=0),
            df2 %>% add_column(join_dummy=0),
            relationship='many-to-many',
            by=c('join_dummy', join_keys),
            ...
        ) %>% 
        select(-c(join_dummy))
    } else {
        full_join(
            df1 %>% add_column(join_dummy=0),
            df1 %>% add_column(join_dummy=0),
            relationship='many-to-many',
            by=c('join_dummy', join_keys),
            ...
        ) %>% 
        select(-c(join_dummy))
    }
}

get_all_row_combinations <- function(
    df1,
    df2=NULL,
    cols_to_pair=c(),
    suffixes=c('.P1', '.P2'),
    keep_self=TRUE,
    ...){
    # Get all combinations of rows with matching attributes (cols_to_pair)
    {
        if (is_tibble(df2)) {
            join_all_rows(
                df1 %>% mutate(index=row_number()),
                df2 %>% mutate(index=row_number()),
                join_keys=cols_to_pair,
                suffix=c('.A', '.B')
            )
        } else {
            join_all_rows(
                df1 %>% mutate(index=row_number()),
                df1 %>% mutate(index=row_number()),
                join_keys=cols_to_pair,
                suffix=c('.A', '.B')
            )
        }
    } %>% 
    # Keep only distinct pairs of rows (order doesnt matter
    mutate(
        index.pair=
            pmap_chr(
                list(index.A, index.B),
                ~ paste(sort(c(...)), collapse='~')
            )
    ) %>%
    distinct(
        index.pair,
        .keep_all=TRUE
    ) %>%
    # keep/remove pairs of a row matched to itself
    {
        if (keep_self) {
            .
        } else {
            filter(., index.A != index.B)
        }
    } %>% 
    # rename with suffixes
    select(-c(starts_with('index'))) %>%
    rename_with(~ str_replace(.x, '\\.A$', suffixes[[1]])) %>% 
    rename_with(~ str_replace(.x, '\\.B$', suffixes[[2]]))
}

merge_sample_info <- function(
    SampleInfo.P1,
    SampleInfo.P2,
    suffix.P1='.P1',
    suffix.P2='.P2',
    ...){
    # SampleInfo.P1=df$SampleInfo.P1[[1]]; prefix.P1='SampleInfo.P1.'; SampleInfo.P2=df$SampleInfo.P2[[1]]; prefix.P2='SampleInfo.P2.'
    # Combine pair info
    if (is.null(SampleInfo.P1)) {
        return(NULL)
    } else if (is.null(SampleInfo.P2)) {
        return(NULL)
    }
    bind_rows(
        SampleInfo.P1 %>% rename_with(~ str_remove(.x, suffix.P1)),
        SampleInfo.P2 %>% rename_with(~ str_remove(.x, suffix.P2)),
    ) %>% 
    add_column(PairIndex=c('P1', 'P2')) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(
        -PairIndex,
        names_to='SampleAttribute',
        values_to='SampleValue'
    ) %>% 
    pivot_wider(
        names_from=PairIndex,
        values_from=SampleValue
    ) %>% 
    rowwise() %>% 
    # Now store both pair values for a single metadata field in 1 column
    # and keep order consistent for grouping downstream
    # NOTE this means that the metadata columns are not always P1 vs P2
    # sometimes it will be P2 vs P1 based on alphabetical order of the values
    # but thats fine since you can just looks at the SampleIDs themselves if 
    # you care about that info
    mutate(
        PairValue=
            case_when(
                P1 == P2 ~ glue('{P1} vs {P2}'),
                P1 != P2 ~ 
                    c(P1, P2) %>% 
                    sort() %>%  
                    paste(collapse=" vs ")
            )
    ) %>%
    select(SampleAttribute, PairValue) %>%
    pivot_wider(
        names_from=SampleAttribute,
        values_from=PairValue
    )
}

###################################################
# Load Specific Data
###################################################
load_sample_metadata <- function(filter=TRUE){
    SAMPLE_METADATA_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    {
        if(filter) {
            filter(., Included)
        } else {
            .
        }
    }
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=c('Chr', 'chr.total.bp')
    )
}

get_min_resolution_per_matrix <- function(
    df,
    as_int=TRUE,
    filter_res=TRUE){
    # get minimum viable resolution for each matrix based on Rao et at. 2014 definition
    MIN_SAMPLE_RESOLUTION_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    select(SampleID, resolution) %>% 
    {
        if (as_int && !is.numeric(.$resolution)) {
            mutate(., resolution=scale_numbers(resolution))
        } else if (!as_int && is.numeric(.$resolution)) {
            mutate(., resolution=scale_numbers(resolution))
        } else if ( as_int &&  is.numeric(.$resolution)) {
            .
        } else if (!as_int && !is.numeric(.$resolution)) {
            .
        }
    } %>% 
    add_column(is.smallest.resolution=TRUE) %>% 
    left_join(
        df,
        by=join_by(SampleID)
    ) %>%
    { 
        if (filter_res) {
            filter(., is.smallest.resolution) %>%
            select(-c(is.smallest.resolution))
        } else {
            .
        }
    }
}

fetch_regions <- function(
    regions.df,
    normalizations,
    resolutions,
    ...){
    regions.df %>% 
    join_all_rows(
        expand_grid(
            normalization=normalizations,
            resolution=resolutions
        )
    ) %>% 
    # add context of 1/2 the size (rounded to nearest bin) to each side of the region
    mutate(
        window.size=
            round(
                region.dist / 2,
                digits=-log10(resolution)
            ),
        across(
            c(
              region.start,
              region.end,
              region.dist,
              resolution,
              window.size
            ),
            as.integer
        )
    ) %>% 
    select(-c(region.group))
}

format_plot_params <- function(
    regions.df,
    title.prefix='RGD Region',
    ...){
    load_mcool_files(
        return_metadata_only=TRUE,
        pattern='*.mapq_30.1000.mcool',
        regions.df=regions.df,
        range1s=NULL,
        range2s=NULL,
        progress=TRUE
    ) %>% 
    # format querys for loading contacts via fetch()
    rowwise() %>% 
    mutate(
        range1=glue('{region.chr}:{max(0, region.start - window.size)}-{region.end + window.size}'),
        range2=range1
    ) %>%
    ungroup() %>% 
    mutate(
        cis=TRUE,
        region.title=glue('{title.prefix} {region} {region.UCSC}'),
        output_dir=
            file.path(
                PLOT_DIR, 
                glue('region_{region}'),
                glue('normalization_{normalization}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('context_{scale_numbers(window.size)}')
            )
    )
}

###################################################
# Load Filetypes
###################################################
load_genome_coverage <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
    rename_with(~ str_remove(.x, '_(raw|weight)')) %>% 
    rename(
        'chr'=chrom,
        'coverage.cis'=cov_cis,
        'coverage.total'=cov_tot,
    ) %>% 
    pivot_longer(
        starts_with('coverage'),
        names_to='metric',
        names_prefix='coverage.',
        values_to='coverage'
    ) %>%  
    mutate(chr=factor(chr, levels=CHROMOSOMES))
}

load_mcool_file <- function(
    filepath,
    resolution=100000,
    normalization="NONE",
    range1="",
    range2="",
    cis=TRUE,
    type='df',
    ...){
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        normalization=normalization,
        join=TRUE,
        query_type='UCSC',
        type=type
    ) %>% 
    as_tibble() %>%
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

list_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    resolutions=NULL,
    normalizations=NULL,
    ...){
    # List all cooler files 
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    # parse sample metadata
    mutate(MatrixID=basename(filepath)) %>% 
    get_info_from_MatrixIDs(
        matrix_ID_col='MatrixID',
        sample_ID_col='SampleID',
        col_prefix='',
        keep_id=FALSE,
        nest_col=NA,
    ) %>% 
    # List all paramter combinations
    {
        if (!is.null(normalizations)) {
            cross_join(., tibble(normalization=normalizations))
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(resolutions)) {
            cross_join(., tibble(resolution=resolutions))
        } else {
            .
        }
    }
}

load_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    resolutions=NULL,
    normalizations=NULL,
    regions.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    return_metadata_only=FALSE,
    keep_metadata_columns=FALSE,
    ...){
    # Define all genomic regions to load contacts 
    regions.df <- 
        if (is.null(regions.df)) {
            # Get the whole genome
            if ((is.null(range1s)) & (is.null(range2s))) {
                tibble()
                # tibble(
                #     range1=CHROMOSOMES,
                #     range2=CHROMOSOMES
                # )
            # get all contacts within all regions in range1s
            } else if (is.null(range2s)) {
                tibble(
                    range1=range1s,
                    range2=range1s
                )
            # get all contacts between all pairs of regions only (not intra-region contacts)
            } else {
                expand_grid(
                    range1=range1s,
                    range2=range2s
                )
            }
        # Just load the specified regions (1 region per row: chr, start, end)
        } else {
            regions.df
        }
    # List all regions for all samples
    list_mcool_files(
        pattern=pattern,
        resolutions=resolutions,
        normalizations=normalizations,
    ) %>% 
    join_all_rows(regions.df) %>% 
    # Load contacts if specified or just return sample metadata + filepaths + regions
    {
        if (return_metadata_only) {
            .
        } else {
            mutate(
                .,
                contacts=
                    purrr::pmap(
                        .l=.,
                        .f=load_mcool_file,
                        .progress=progress
                    )
            ) %>% 
            select(-c(filepath, range1, range2)) %>% 
            unnest(contacts)
        }
    }
}

