library(tidyverse)
library(magrittr)
library(tictoc)
library(glue)
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
    mutate(Sample.ID=glue('{Edit}.{Genotype}.{SampleNumber}.{Celltype}')) %>% 
    select(-c(Edit, Genotype, SampleNumber, Celltype))
}

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
    ...){
    filepath %>% 
    CoolFile() %>% 
    import(...)
}

load_sparse_matrix <- function(filename){
    read.table(
        filename,
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

load_all_sparse_matrix_files <- function(
    matrix_dir,
    ...){
    matrix_dir %>%
    list.files(
        pattern='*.txt',
        full.names=FALSE
    ) %>%
    tibble(filename=.) %>% 
    mutate(
        filepath=file.path(matrix_dir, filename),
        filename=
            filename %>% 
            str_remove('.txt')
    ) %>% 
    separate_wider_delim(
        filename,
        delim=fixed('.'),
        names=
            c(
                'SampleID',
                NA,
                'ReadFilter',
                NA,
                'Resolution',
                'Chromosome'
            )
    ) 
}

fetch_contacts <- function(
    mcool_file,
    resolution,
    normalization,
    ...){
    # hic.obj=hic_matrices$matrix[[1]]; SampleID=hic_matrices$SampleID[[1]]
    interactions.file <- mcool_file %>% File(resolution=resolution)
    cnv.interactions <- 
        interactions.file %>%
        fetch(
            range1=CHR16_TELOMERE_REGION,
            range2=CHR16_ELISE_REGION,
            normalization=normalization,
            join=TRUE,
            query_type='UCSC',
            type='df'
        ) %>% 
        as_tibble() %>% 
        add_column(region='16p11.2 CNV - Telomere')
    dist.range <- 
        cnv.interactions %>%
        mutate(dist=start2 - start1) %>%
        {list('min'=min(.$dist), 'max'=max(.$dist))}
    telomere <- 
        CHR16_TELOMERE_REGION %>% 
        str_remove('^chr16:') %>% 
        str_remove_all(',') %>% 
        str_split('-') %>% 
        first() 
    # Get Telomere+CNV boundaries to not double-count them in controls
    telomere.start <- as.integer(telomere[[1]])
    telomere.end <- as.integer(telomere[[2]])
    cnv <- 
        CHR16_ELISE_REGION %>% 
        str_remove('chr16:') %>% 
        str_remove_all(',') %>% 
        str_split('-') %>% 
        first()
    cnv.start <- as.integer(cnv[[1]])
    cnv.end <- as.integer(cnv[[2]])
    # Get sets of control interactions to compare against i.e.
    # all contacts within and distance range as Telomere-CNV contacts
    CHR16_CONTROL_REGIONS %>% 
    as_tibble() %>%
    # First get contacts within each region of interest
    pivot_longer(
        everything(),
        names_to='region',
        values_to='coords'
    ) %>%
    # Now find all contacts within specified distance range
    mutate(
        region=glue('{region}.Distance.Matched.Control'),
        interactions=
            pmap(
                .l=.,
                .f=function(interactions.file, coords, ...){
                    interactions.file %>%
                    fetch(
                        range1=coords,
                        normalization=normalization,
                        join=TRUE,
                        query_type='UCSC',
                        type='df'
                    )
                },
                interactions.file=interactions.file
            )
    ) %>% 
    unnest(interactions) %>% 
    mutate(dist=start2 - start1) %>%
    dplyr::filter(between(dist, dist.range$min, dist.range$max)) %>%
    # Remove any contacts spanning regions of interest
    filter(
        !(
            between(start1, telomere.start, telomere.end) &
            between(start2, cnv.start, cnv.end)
        )
    ) %>% 
    select(-dist) %>% 
    # add cnv contacts back and return 
    bind_rows(cnv.interactions) %>%
    mutate(
        region=
            factor(
                region,
                levels=
                    CHR16_CONTROL_REGIONS %>% 
                    names() %>% 
                    paste0('.Distance.Matched.Control') %>% 
                    sort() %>% 
                    c('16p11.2 CNV - Telomere', .)
            )
    ) %>% 
    select(-c(coords, end1, end2))
}
