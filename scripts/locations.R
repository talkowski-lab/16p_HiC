library(here)
library(tidyverse)
library(magrittr)
library(glue)
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
# COVERAGE_DATA_FILE <-    file.path(COVERAGE_DIR, 'all.coverage.data.tsv')
# hicrep results
HICREP_DIR <-            file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE <-   file.path(HICREP_DIR, 'all.hicrep.scores.tsv')
# Reproducing figure from Elise Robinson
ELISE_RECREATION_DIR <-  file.path(RESULTS_DIR, 'Elise.Recreation.results')
# Differential Contact results from multiHiCCompare
SPARSE_MATRIX_DIR <-     file.path(RESULTS_DIR, 'sparse.matrices')
MULTIHICCOMPARE_DIR <-   file.path(RESULTS_DIR, 'multiHiCCompare')
# TAD Annotations
TAD_DIR <-               file.path(RESULTS_DIR, 'TADs')
TAD_RESULTS_FILE <-      file.path(TAD_DIR, 'all.TAD.annotations.tsv')
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
    tibble(resolution=.) %>% 
    mutate(resolution.name=resolution / 1000) %>% 
    mutate(
        resolution.name=
            case_when(
                resolution.name >= 1000 ~ paste0(resolution.name / 1000, "Mb"),
                resolution.name >= 1    ~ paste0(resolution.name, "Kb")
            ) %>% 
            factor(., levels=.)
    )
# Match ideal smoothing param to specified resolution based on this
# https://github.com/TaoYang-dev/hicrep?tab=readme-ov-file#hicrep-parameters
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
# TO have consistent colors for genotypes across figures
GENOTYPE_COLORS <- 
    c(
        'WT'='#b3b3ff',
        'DEL'='#ff0000',
        'DUP'='#0000ff'
    )
# Genomic region boundaries to refernce during analysis plotting
GENOMIC_REGIONS <- 
    tribble(
        ~region,           ~region.chr, ~region.start, ~region.end,
        "chr16p.deletion",     "chr16",      29488679,    30188679,
        "chr16p.telomere",     "chr16",             0,     5149999,
        "chr16p11.2",          "chr16",      24288679,    30188679,
        "chr16p",              "chr16",             0,    36800000,
        "chr16",               "chr16",             0,    90338345,
        "NIPBL",                "chr5",      36876769,    37066413,
        "chr5p13.2",            "chr5",      33800001,    38400000,
        "chr5p",                "chr5",             0,    48800000,
        "chr5",                 "chr5",             0,   181538259,
        "WAPL",                "chr10",      86435256,    86521792,
        "chr10q23.2",          "chr10",      86100001,    87700000,
        "chr10q",              "chr10",      39800001,   133797422,
        "chr10",               "chr10",             0,   133797422,
    ) %>% 
    mutate(region.UCSC=glue("{region.chr}:{region.start}-{region.end}")) %>%
    mutate(region.dist=region.end - region.start)
