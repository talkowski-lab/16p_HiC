library(magrittr)
library(tictoc)
library(tidyverse)
###############
# Load Data
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
