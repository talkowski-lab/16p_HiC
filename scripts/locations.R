library(here)
library(tidyverse)
library(magrittr)
BASE_DIR <- here()
# Input Data
GENOME_REF_DIR <-       '/data/talkowski/tools/ref/Hi_c_noalt'
CHROMOSOME_SIZES_FILE <- file.path(GENOME_REF_DIR, 'hg38.reduced.chrom.sizes')
FASTQ_DIR <-            file.path(BASE_DIR, 'fastq')
SAMPLE_METADATA_FILE <- file.path(BASE_DIR, 'sample.metadata.tsv')
# Results Dirs
RESULTS_DIR <-          file.path(BASE_DIR, 'results')
PAIRS_DIR <-            file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
COOLERS_DIR <-          file.path(RESULTS_DIR, 'coolers_library')
HICREP_DIR <-           file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE <-  file.path(HICREP_DIR, 'all.hicrep.scores.tsv')
MATRIX_QC_DIR <-        file.path(RESULTS_DIR, 'sample.QC')
ELISE_RECREATION_DIR <- file.path(RESULTS_DIR, 'Elise.Recreation.results')
TAD_DIR <-              file.path(RESULTS_DIR, 'TADs')
# SPARSE_MATRIX_DIR <-    file.path(RESULTS_DIR, 'sparse.matrices')
# MULTIHICCOMPARE_DIR <-  file.path(RESULTS_DIR, 'multiHiCCompare')
# Factor levels for variaous metadata categories
EDITS <- 
    c(
      '16p',
      'WAPL',
      'NIBPL',
      'RAD21'
    )
GENOTYPES <- 
    c(
      'WT',
      'DEL',
      'DUP'
    )
CELLTYPES <- 
    c(
      'NSC',
      'iN'
    )
CHROMOSOMES <- 
    c(
        paste0('chr', 1:22), 
        'chrX',
        'chrY'
    )
RESOLUTION_NAMES <- 
    tribble(
       ~Resolution, ~Resolution.name, 
              5000,            '5Kb',
             10000,           '10Kb',
             25000,           '25Kb',
             40000,           '40Kb',
             50000,           '50Kb',
            100000,          '100Kb',
            500000,          '500Kb',
           1000000,            '1Mb'
    ) %>%
    mutate(Resolution.name=factor(Resolution.name, levels=Resolution.name))
RESOLUTIONS <- 
    RESOLUTION_NAMES$RESOLUTIONS
RESOLUTION_IDEAL_H <- 
    tribble(
       ~Resolution, ~Ideal_H, 
             10000,       20,
             25000,       10,
             40000,        5,
            100000,        3,
            500000,        2,
            500000,        1,
           1000000,        1,
           1000000,        0
    )
