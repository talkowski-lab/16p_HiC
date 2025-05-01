library(here)
library(tidyverse)
library(magrittr)
BASE_DIR <- here()
# Input Data
GENOME_REF_DIR <-       '/data/talkowski/tools/ref/Hi_c_noalt'
CHROMOSOME_SIZES_FILE <- file.path(GENOME_REF_DIR, 'hg38.reduced.chrom.sizes')
FASTQ_DIR <-             file.path(BASE_DIR, 'fastq')
SAMPLE_METADATA_FILE <-  file.path(BASE_DIR, 'sample.metadata.tsv')
# distiller produced ou  puts
RESULTS_DIR <-           file.path(BASE_DIR, 'results')
BAM_DIR <-               file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
PAIRS_DIR <-             file.path(RESULTS_DIR, 'pairs_library')
COOLERS_DIR <-           file.path(RESULTS_DIR, 'coolers_library')
# Used for Sample QCing  
SAMPLE_QC_DIR <-         file.path(RESULTS_DIR, 'sample.QC')
COVERAGE_DIR <-          file.path(SAMPLE_QC_DIR, 'coverage')
COVERAGE_DATA_FILE <-    file.path(COVERAGE_DIR, 'all.coverage.data.tsv')
# hicrep results
HICREP_DIR <-            file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE <-   file.path(HICREP_DIR, 'all.hicrep.scores.tsv')
# Reproducing figure fr  m Elise Robinson
ELISE_RECREATION_DIR <-  file.path(RESULTS_DIR, 'Elise.Recreation.results')
# Differential Contact   esults from multiHiCCompare
MULTIHICCOMPARE_DIR <-   file.path(RESULTS_DIR, 'multiHiCCompare')
SPARSE_MATRIX_DIR <-     file.path(MULTIHICCOMPARE_DIR, 'sparse.matrices')
# TAD Annotations
TAD_DIR <-               file.path(RESULTS_DIR, 'TADs')
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
RESOLUTIONS <- 
    c(
          1000,
          5000,
         10000,
         25000,
         40000,
         50000,
        100000,
        500000,
       1000000,
       2500000,
       5000000
    )
RESOLUTION_NAMES <- 
    RESOLUTIONS %>%
    tibble(Resolution=.) %>% 
    mutate(Resolution.name=Resolution / 1000) %>% 
    mutate(
        Resolution.name=
            case_when(
                Resolution.name >= 1000 ~ paste0(Resolution.name / 1000, "Mb"),
                Resolution.name >= 1    ~ paste0(Resolution.name, "Kb")
            ) %>% 
            factor(., levels=.)
    )
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
GENOTYPE_COLORS <- 
    c(
        'WT'='#b3b3ff',
        'DEL'='#ff0000',
        'DUP'='#0000ff'
    )
# Genomic region boundaries to refernce during analysis plotting
GENOMIC_REGIONS <- 
    tribble(
        ~name,              ~chr,    ~start,       ~end,
        "chr16",              16,         0,   90338345,
        "chr16p",             16,         0,   36800000,
        "chr16p11.2",         16,  24288679,   30188679,
        "chr16p.deletion",    16,  29488679,   30188679,
        "chr16p.telomere",    16,         0,    5149999,
        "chr5",                5,         0,  181538259,
        "chr5p",               5,         0,   48800000,
        "chr5p13.2",           5,  33800001,   38400000,
        "NIPBL",               5,  36876769,   37066413,
        "chr10",              10,         0,  133797422,
        "chr10q",             10,  39800001,  133797422,
        "chr10q23.2",         10,  86100001,   87700000,
        "WAPL",               10,  86435256,   86521792
    ) %>% 
    mutate(UCSC=glue("{chr}:{start}-{end}"))
