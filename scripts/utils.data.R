library(tidyverse)
library(magrittr)
library(tictoc)
library(glue)
library(HiCExperiment)
library(hictkR)
###############
# Pairsing and Caching
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
    keep_metadata_columns=FALSE,
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
    {
        if (keep_metadata_columns) {
            .
        } else {
            select(., -c(Edit, Genotype, SampleNumber, Celltype))
        }
    }
}

###############
# Utilities
scale_numbers <- function(
    numbers,
    accuracy=2){
    print(is.numeric(numbers))
    if (is.character(numbers)) {
        case_when(
            grepl('Mb', numbers) ~ as.integer(str_remove(numbers, 'Mb')) * 1e6,
            grepl('Kb', numbers) ~ as.integer(str_remove(numbers, 'Kb')) * 1e3,
            TRUE ~ as.integer(numbers)
        )
    } else if (is.numeric(numbers)) {
        case_when(
            numbers >= 1e6 ~ glue('{signif(numbers / 1e6, digits=accuracy)}Mb'),
            numbers >= 1e3 ~ glue('{signif(numbers / 1e3, digits=accuracy)}Kb'),
            TRUE ~ as.character(numbers)
        )
    } else {
        numbers
    }
}

count_dirs <- function(filepath){
    filepath %>%
    dirname() %>% 
    str_split('/') %>% 
    length()
}

join_all_rows <- function(
    df1,
    df2,
    ...){
    full_join(
        df1 %>% add_column(join_dummy=0),
        df2 %>% add_column(join_dummy=0),
        relationship='many-to-many',
        by=join_by(join_dummy),
        ...
    ) %>% 
    select(-c(join_dummy))
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
        # all possible pairs of regions to annotate
        # each bin-pair will be assigned 1 region pair, even if regions overlap/contain eachother
        join_all_rows(
            .,
            {.},
            suffix=c('1','2'),
            relationship='many-to-many',
        ) %>% 
        filter(region.chr1 == region.chr2) %>% 
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

fetch_RGD_regions <- function(
    normalizations,
    resolutions,
    ...){
    GENOMIC_REGIONS %>% 
    filter(region.group == "RGDs") %>% 
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

format_rgd_plot_params <- function(
    region.df,
    ...){
    load_mcool_files(
        return_metadata_only=TRUE,
        pattern='*.mapq_30.1000.mcool',
        region.df=region.df,
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
        region.title=glue('{region} RGD Region {region.UCSC}'),
        output_dir=
            file.path(
                PLOT_DIR, 
                glue('region_{region}'),
                glue('normalization_{normalization}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('context_{scale_numbers(window.size)}')
            )
    )
    # group_by(Sample.ID, region) %>% 
    # add_tally(wt=IF, name='region.total.contacts') %>% 
    # ungroup() %>%
    # group_by(Sample.ID, region) %>% 
    # add_count(name='region.nbins') %>% 
    # ungroup() %>%
}
###############
# Load Filetypes
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
    {
        if (normalization == 'weight') {
            add_column(., weight='ICE')
        } else {
            add_column(., weight=normalization)
        }
    } %>% 
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
    region.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    return_metadata_only=FALSE,
    keep_metadata_columns=FALSE,
    ...){
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    mutate(matrix.name=basename(filepath)) %>% 
    {
        if (!is.null(region.df)) {
            join_all_rows(
                .,
                region.df
            )
        } else {
            if ((is.null(range1s)) & (is.null(range2s))) {
                .
            } else if (is.null(range2s)) {
                join_all_rows(
                    .,
                    tibble(
                        range1=range1s,
                        range2=range1s
                    )
                )
            } else {
                join_all_rows(
                    .,
                    expand_grid(
                        range1=range1s,
                        range2=range2s
                    )
                )
            }
        }
    } %>% 
    process_matrix_name(keep_metadata_columns=keep_metadata_columns) %>% 
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
