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
    # This shouldnt  break if the filename has a param_delim character in it 
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

load_mcool_file <- function(
    filepath,
    resolution,
    range1="",
    range2="",
    normalization='NONE',
    cis=TRUE,
    ...){
    # filepath %>% 
    # CoolFile(resolution=resolution) %>%
    # import(focus=region)
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        join=TRUE,
        query_type='UCSC',
        type='df'
    ) %>% 
    as_tibble() %>%
    add_column(weight=normalization) %>% 
    # cis only
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
    range1s,
    range2s=NULL,
    ...){
    # resolutions=c('100000'); regions='XVI'; pattern='mapq_30.1000.mcool' 
    range2s <- ifelse(is.null(range2s), range1s, range2s)
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    mutate(matrix.name=basename(filepath)) %>% 
    mutate(join_dummy=0) %>% 
    full_join(
        expand_grid(
            resolution=resolutions,
            range1=range1s,
            range2=range2s,
            join_dummy=0
        )
    ) %>% 
    process_matrix_name() %>% 
    mutate(
        contacts=
            purrr::pmap(
                .l=.,
                .f=load_mcool_file,
                .progress=TRUE
            )
    ) %>% 
    mutate(resolution=as.integer(resolution)) %>% 
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
