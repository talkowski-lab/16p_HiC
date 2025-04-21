library(magrittr)
library(tictoc)
library(tidyverse)
###############
# Misc
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    results_fnc,
    ...){
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
    # Now check if the results file exists and load it
    tic()
    if (file.exists(results_file) & !force_redo) {
        message(glue('{results_file} exists, not recomputing results'))
        if (return_data) {
            message('Loading cached results...')
            results <- load_fnc(results_file)
        }
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
    toc()
    if (return_data) return(results) else return(invisible(NULL))
}
###############
# Load Data
load_sample_metadata <- function(clean=TRUE){
    SAMPLE_METADATA_FILE %>%
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
# Matrix Utils
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

summarize_contacts_per_bin <- function(
    mcool_file,
    resolution,
    normalization,
    cluster=NULL,
    ...){
    # All pairs of bins with > 0 raw count as df
    interactions.df <- 
        mcool_file %>% 
        File(resolution=resolution) %>%
        fetch(
            normalization=normalization,
            join=TRUE,
            query_type='UCSD',
            type='df'
        ) %>% 
        as_tibble() %>% 
        select(
            chrom1,
            start1,
            chrom2,
            start2,
            count
        )
    # bin_x == Left side
    A.interactions <- 
        interactions.df %>%
        rename(
            chrom1='id.chr',
            start1='id.bin',
            chrom2='nest.chr',
            start2='nest.bin'
        ) %>% 
        group_by(id.chr, id.bin) %>%
        nest(contacts=c(nest.chr, nest.bin, count)) %>%
        ungroup()
    # bin_x == Rightside
    B.interactions <- 
        interactions.df %>%
        rename(
            chrom2='id.chr',
            start2='id.bin',
            chrom1='nest.chr',
            start1='nest.bin'
        ) %>% 
        group_by(id.chr, id.bin) %>%
        nest(contacts=c(nest.chr, nest.bin, count)) %>%
        ungroup()
    #  Now, to be sure each bin e.g. bin_x is tallied correctly, do this
    #  1. get frequency table for all contacts where left  bin == bin_x
    #  2. get frequency table for all contacts where right bin == bin_x
    #  3. concat both tables for bin_x i.e. # of times bin_x forms a contact
    inner_join(
        A.interactions, B.interactions, by=join_by(id.chr, id.bin),
        suffix=c('.A', '.B')
    ) %>%
    # Compute summary stats per bin for all observed pairs
    rowwise() %>% 
    mutate(
        contacts=
            bind_rows(
                contacts.A,
                contacts.B
            ) %>% 
            list()
    ) %>% 
    select(-c(contacts.A, contacts.B)) %>% 
    unnest(contacts) %>%
    group_by(id.chr, id.bin) %>% 
    { 
        if (is.null(cluster)) {
            .
        } else {
            partition(., cluster)
        }
    } %>% 
    summarize(
        n=n(),
        var=var(count, na.rm=TRUE),
        total=sum(count, na.rm=TRUE),
        summary=list(as_tibble_row(summary(count, na.rm=TRUE)))
    ) %>% 
    { 
        if (is.null(cluster)) {
            .
        } else {
            collect(.)
        }
    } %>% 
    ungroup() %>%
    unnest(summary) %>% 
    mutate(across(where(~ class(.x) == 'table'), ~ as.numeric(.x))) %>% 
    # clean before saving
    pivot_longer(
        -c(id.chr, id.bin),
        names_to='stat',
        values_to='value',
    ) %>%
    rename(id.chr='Chr', id.bin='bin.start') %>% 
    arrange(Chr, bin.start)
}
load_pairtools_stats <- function(
    pairs_dir,
    ...){
    # Load all stats
    all.stats.df <- 
        pairs_dir %>%
        list.files(
            full.names=TRUE,
            # pattern='*.dedup.stats'
            pattern='*.pairs_stats.tsv'
        ) %>%
        tibble(filepath=.) %>% 
        mutate(filename=basename(filepath)) %>% 
        separate_wider_delim(
            filename,
            delim=fixed('.'),
            names=
                c(
                    'SampleID',
                    NA,
                    NA,
                    NA
                )
        ) %>%
        mutate(
            stats=
                purrr::pmap(
                    .l=.,
                    function(filepath, ...){
                        filepath %>% 
                        read_tsv(
                            show_col_types=FALSE,
                            col_names=c('stat', 'value')
                        )
                    }
                )
        ) %>%
        select(-c(filepath)) %>%
        unnest(stats)
    # Matrix-wide stats
    general.stats.df <- 
        all.stats.df %>% 
        filter(!grepl('/', stat)) %>%
        mutate(
            category=
                case_when(
                    grepl('total', stat) ~ 'Reads', 
                    TRUE ~ 'Contacts'
                )
        ) 
    # Get total number of contacts per sample
    total.contacts.df <- 
        general.stats.df %>% 
        filter(stat == 'total_nodups') %>% 
        select(SampleID, value) %>%
        dplyr::rename('total.unique.contacts'=value)
    general.stats.df <- 
        general.stats.df %>% left_join(total.contacts.df, by='SampleID')
    # stats about contact distances
    dist.stats.df <- 
        all.stats.df %>% 
        filter(grepl('^dist_freq', stat)) %>%
        separate_wider_delim(
            stat,
            delim='/',
            names=
                c(
                    'interaction',
                    'range',
                    'orientation'
                )
        ) %>%
        mutate(range=str_remove(range, fixed('+'))) %>% 
        separate_wider_delim(
            range,
            delim='-',
            names=c('range.start', 'range.end'),
            too_few='align_start'
        ) %>%
        mutate(across(starts_with('range'), ~ as.integer(.x))) %>% 
        mutate(interaction='cis') %>%
        left_join(total.contacts.df, by='SampleID')
    # Contact Frequency Stats per Chromosome
    chr.stats.df <- 
        all.stats.df %>% 
        filter(grepl('^chrom_freq', stat)) %>%
        separate_wider_delim(
            stat,
            delim='/',
            names=
                c(
                    'stat',
                    'chr1',
                    'chr2'
                )
        ) %>% 
        left_join(total.contacts.df, by='SampleID')
    # Return named list of all stuff
    return(
        list(
            'general.stats'=general.stats.df,
            'distance.stats'=dist.stats.df,
            'chr.stats'=chr.stats.df
        )
    )
}
