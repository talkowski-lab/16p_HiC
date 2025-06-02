library(here)
library(tidyverse)
library(magrittr)
library(glue)
BASE_DIR <- here()
# Input Data
GENOME_REF_DIR                 <- '/data/talkowski/tools/ref/Hi_c_noalt'
CHROMOSOME_SIZES_FILE          <- file.path(GENOME_REF_DIR, 'hg38.reduced.chrom.sizes')
FASTQ_DIR                      <- file.path(BASE_DIR, 'fastq')
SAMPLE_METADATA_FILE           <- file.path(BASE_DIR, 'sample.metadata.tsv')
# distiller produced ou           puts
RESULTS_DIR                    <- file.path(BASE_DIR, 'results')
BAM_DIR                        <- file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
PAIRS_DIR                      <- file.path(RESULTS_DIR, 'pairs_library')
COOLERS_DIR                    <- file.path(RESULTS_DIR, 'coolers_library')
# Used for Sample QCing           
SAMPLE_QC_DIR                  <- file.path(RESULTS_DIR, 'sample.QC')
COVERAGE_DIR                   <- file.path(SAMPLE_QC_DIR, 'coverage')
# COVERAGE_DATA_FILE <-             file.path(COVERAGE_DIR, 'all.coverage.data.tsv')
# hicrep results
HICREP_DIR                     <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE            <- file.path(HICREP_DIR, 'all.hicrep.scores.tsv')
# Reproducing figure from Elise Robinson
ROBINSON_REPLICATION_DIR       <- file.path(RESULTS_DIR, 'Robinson.replication')
ROBINSON_REPLICATION_DATA_FILE <- file.path(ROBINSON_REPLICATION_DIR, 'replication.data.chr16.tsv')
# Differential Contact results from multiHiCCompare
SPARSE_MATRIX_DIR              <- file.path(RESULTS_DIR, 'sparse.matrices')
MULTIHICCOMPARE_DIR            <- file.path(RESULTS_DIR, 'multiHiCCompare')
# TAD Annotations
TAD_DIR                        <- file.path(RESULTS_DIR, 'TADs')
TAD_RESULTS_FILE               <- file.path(TAD_DIR, 'all.TAD.annotations.tsv')
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
# 16p deleteion and telomere are to replicate Fig 4B in this paper:
# https://www.nature.com/articles/s41588-022-01203-y
GENOMIC_REGIONS <- 
	tribble(
    ~region.group, ~region,                    ~region.chr, ~region.start, ~region.end,
    "16p",         "chr16p.deletion",          "chr16",          29488679,    30188679,
    "16p",         "chr16p.telomere",          "chr16",                 0,     5149999,
    "16p",         "chr16p11.2",               "chr16",          24288679,    30188679,
    "16p",         "chr16p",                   "chr16",                 0,    36800000,
    "16p",         "chr16q",                   "chr16",          36800001,    90338345,
    "16p",         "chr16",                    "chr16",                 0,    90338345,
    "NIPBL",       "NIPBL",                    "chr5",           36876769,    37066413,
    "NIPBL",       "chr5p13.2",                "chr5",           33800001,    38400000,
    "NIPBL",       "chr5p",                    "chr5",                  0,    48800000,
    "NIPBL",       "chr5",                     "chr5",                  0,    81538259,
    "WAPL",        "WAPL",                     "chr10",          86435256,    86521792,
    "WAPL",        "chr10q23.2",               "chr10",          86100001,    87700000,
    "WAPL",        "chr10q",                   "chr10",          39800001,    33797422,
    "WAPL",        "chr10",                    "chr10",                 0,    33797422,
    "RGDs",        "1q21",                     "chr1",          146081967,   148779515,
    "RGDs",        "7q11.23_WBS",              "chr7",           73174795,    74833380,
    "RGDs",        "8p23.1",                   "chr8",             523442,     5539949,
    "RGDs",        "15q_AS_PWS_large",         "chr15",          22420444,    28730223,
    "RGDs",        "15q13.3",                  "chr15",          30377807,    32441361,
    "RGDs",        "16p13.11",                 "chr16",          14816373,    16353400,
    "RGDs",        "16p11.2",                  "chr16",          29449194,    30335547,
    "RGDs",        "17q12",                    "chr17",          36222499,    38194356,
    "RGDs",        "17q21.31",                 "chr17",            452200,     1272274,
    "RGDs",        "22q11.21_DGS_VCFS_common", "chr22",          18518837,    21562827,
	)	%>%	
    mutate(
        # region=fct_reorder(region, region.dist)
        region.UCSC=glue("{region.chr}:{region.start}-{region.end}"),
        region.dist=region.end - region.start
    )
