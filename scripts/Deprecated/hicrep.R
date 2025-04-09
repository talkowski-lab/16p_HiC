library(glue)
library(tidyverse)
library(furrr)
library(hicrep)
# Dirs
# BASE_DIR='/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p'
BASE_DIR='/data/talkowski/broadIncoming/22LCC2LT4/fastq'
source(glue('{BASE_DIR}/scripts/locations.R'))
source(glue('{BASE_DIR}/scripts/utils.data.R'))
# Matrix data
SPARSE_DIR=glue('{BASE_DIR}/results/sparse.matrices')
# COOLER_DIR=glue('{BASE_DIR}/results/coolers_library')
CHROMOSOMES=c(paste('chr', 1:22), 'chrX', 'chrY')
# list matrix files
matrix_files <- 
    SPARSE_DIR %>%
    load_all_matrix_files() %>% 
    mutate(
        Resolution=as.integer(Resolution),
        Genotype=
            SampleID %>% 
            str_extract('16p(DEL|DUP|WT)') %>%
            str_remove('16p')
    )
matrix_files %>% dplyr::count(SampleID)
matrix_files %>% dplyr::count(Resolution, Genotype, ReadFilter, SampleID)
matrix_files %>% select(-c(filepath)) %>% head(3)
# From hicrep github https://github.com/TaoYang-dev/hicrep
h_ideals <- 
    tribble(
        ~Resolution, ~h_ideal,
       10000,20,
       25000,10,
       40000,5,
       100000,3,
       500000,2,
       1000000,0
    )
# Run hicrep 
hicrep_results <- 
    # generate all pairs of samples matched in relevant variables
    full_join(
        x=matrix_files,
        y=matrix_files,
        by=join_by(Chromosome, Resolution, ReadFilter),
        suffix=c('.1', '.2'),
        relationship='many-to-many',
        multiple='all'
    ) %>%
    left_join(
        h_ideals,
        by=join_by(Resolution),
        multiple='all'
    )
    hicrep_results %>% 
    filter(Resolution == 10000) %>%
    filter(Chromosome == 'chr1') %>%
    # select(-c(filepath.1, filepath.2)) %>% head(3) %>% t()
    # dplyr::count(ReadFilter, SampleID.1)
    # future_pmap(
    head(2) %>% tail(1) %>% 
    pmap(
        .l=.,
        .f=run_hicrep,
        .progress=TRUE
    )
DEPRECIATED <- function(){
    h_ideals <- 
        matrix_files %>%
        group_by(Resolution, ReadFilter) %>% 
        filter(Chromosome == 'chr10') %>% 
        mutate(h_train_pair=glue('{SampleID}:{Chromosome}')) %>% 
        sample_n(2) %>% 
        mutate(tmp=row_number()) %>% 
        pivot_wider(
            id_cols=c(Resolution, ReadFilter),
            names_prefix='filepath.',
            names_from=tmp,
            values_from=filepath
        ) %>%
        ungroup() %>% 
        # head(1) %>% t()
        head(2) %>% future_pmap(.l=., .f=run_htrain)
        mutate(
            h_ideal=
                future_pmap(
                    .l=.,
                    .f=run_htrain,
                )
        )
    }
