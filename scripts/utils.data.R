suppressMessages(library(magrittr))
suppressMessages(library(tictoc))
suppressMessages(library(tidyverse))
source(glue('{BASE_DIR}/scripts/locations.R'))
source(glue('{BASE_DIR}/scripts/utils.plot.R'))
source(glue('{BASE_DIR}/scripts/utils.TADs.R'))
CHROMOSOMES <- c(paste0('chr', 1:22), 'chrX', 'chrY', 'chrM')
RESOLUTION_NAMES <- 
    tribble(
       ~ Resolution, ~ Resolution.name, 
       '5000','5Kb',
       '10000','10Kb',
       '50000','50Kb',
       '100000','100Kb',
       '500000','500Kb',
       '1000000','1Mb'
    ) %>%
    mutate(Resolution.name=factor(Resolution.name, c('5Kb', '10Kb', '50Kb', '100Kb', '500Kb', '1Mb')))
       
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
# Sample Metadata
get_genotype <- function(SampleID){
    case_when(
        str_detect(SampleID, '^16pWT')  ~ 'WT',
        str_detect(SampleID, '^16pDEL') ~ 'DEL',
        str_detect(SampleID, '^16pDUP') ~ 'DUP',
        TRUE ~ 'Unknown'
    )
}
get_celltype <- function(SampleID){
    case_when(
        str_detect(SampleID, 'NSCHiC$') ~ 'NSC',
        str_detect(SampleID, 'iNHiC$')   ~ 'iN',
        TRUE ~ 'Unknown'
    )
}
make_sample_metadata <- function(output_file=SAMPLE_METADATA_TSV){
    FASTQ_DIR  %>% 
    # List all fastq files (2 per sample)
    list.files(pattern='*16p*.fastq.gz') %>%
    tibble(filename=.) %>%
    mutate(fastq=glue('{FASTQ_DIR}/{filename}')) %>%
    # extact info from filename
    separate_wider_delim(
        filename,
        delim='_',
        too_few='align_start',
        names=
            c(
                'Sequencing.Run',
                'Lane',
                'Investigator.ID',
                'Sample.Name',
                NA,
                NA,
                'fastq_side',
                NA
            )
    ) %>%
    # 1 row per sample, with R1,R2 files as columns
    pivot_wider(
        names_from=fastq_side,
        values_from=fastq 
    ) %>% 
    # Get date experiments were submitted
    left_join(
        glue('{BASE_DIR}/dragen.stats/sample.metadata.tsv') %>%
        read_csv() %>%
        select(c(Sample_ID, Sample_Name, Date)),
        by=join_by(Sample.Name == Sample_Name)
    ) %>% 
    rename('Sample.ID'=Sample_ID) %>%
    # extract more info from Sample Name
    mutate(
        CellType=
            get_celltype(Sample.Name) %>% 
            factor(levels=c('NSC', 'iN')),
        Genotype=
            get_genotype(Sample.Name) %>% 
            factor(levels=c('WT', 'DEL', 'DUP'))
    ) %>%
    write_tsv(output_file)
}
load_sample_metadata <- function(clean=TRUE){
    COOLERS_DIR %>%
    list.files(
        pattern='*.mapq_30.1000.mcool',
        full.names=TRUE
    ) %>% 
    tibble(mcool_file=.) %>%
    mutate(
        filename=basename(mcool_file),
        Sample.Matrix=str_remove(filename, '.mcool')
    ) %>% 
    separate_wider_delim(
        filename,
        delim=fixed('.'),
        names=c('Sample.ID', NA, 'ReadFilter', NA, NA)
    ) %>%
    mutate(
        is.merged=grepl('S[0-9]S[0-9]$', Sample.ID),
        Replicate.ID=
            case_when(
                is.merged ~ 'Merged',
                TRUE ~ Sample.ID %>%
                       str_extract('^16p.*_S([1-6])\\.*', group=1) %>% 
                       as.integer() %>% 
                       { 2 - . %% 2 } %>%
                       as.character() %>%
                       paste0('Rep', .)
            ) %>% 
            factor(levels=c('Rep1', 'Rep2', 'Merged')),
        Sample.ID=str_remove(Sample.ID, '_S[0-9].*$')
    ) %>% 
    left_join(
        read_tsv(SAMPLE_METADATA_TSV),
        by=join_by(Sample.ID == Sample.Name)
    ) %>%
    {
        if (clean) {
            select(., -c(R1, R2, Lane, Investigator.ID))
        } else {
            .
        }
    }
}
###############
# Load Data
load_chr_sizes <- function(){
    '/data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes' %>%
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
        filepath=glue('{matrix_dir}/{filename}'),
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
    ) #%>%
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
        A.interactions,
        B.interactions,
        by=join_by(id.chr, id.bin),
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
###############
# DEPRECATED
DEP_load_16p_sample_metadata <- function(){
    COOLERS_DIR %>%
    list.files(
        pattern='*.mapq_30.1000.mcool',
        full.names=TRUE
    ) %>% 
    tibble(mcool_file=.) %>%
    mutate(
        filename=basename(mcool_file),
        SampleMatrix=str_remove(filename, '.mcool')
    ) %>% 
    separate_wider_delim(
        filename,
        delim=fixed('.'),
        names=c('SampleID', NA, 'ReadFilter', NA, NA)
    ) %>%
    mutate(
        Genotype=
            case_when(
                str_detect(SampleID, '^16pWT') ~ 'WT',
                str_detect(SampleID, '^16pDEL') ~ 'DEL',
                str_detect(SampleID, '^16pDUP') ~ 'DUP',
                TRUE ~ 'Unknown'
            ) %>%
            factor(levels=c('WT', 'DEL', 'DUP')),
        is.merged=grepl('S[0-9]S[0-9]$', SampleID),
        Replicate.ID=
            case_when(
                is.merged ~ 'Merged',
                TRUE ~ SampleID %>%
                       str_extract('^16p.*_S([1-6])\\.*', group=1) %>% 
                       as.integer() %>% 
                       { 2 - . %% 2 } %>%
                       as.character() %>%
                       paste0('Rep', .)
            ) %>% 
            factor(levels=c('Rep1', 'Rep2', 'Merged'))
    )
}
