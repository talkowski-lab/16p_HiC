# Dependencies
library(tidyverse)
library(magrittr)
library(glue)
# Factor levels for variaous metadata categories
CNV.GROUPS <- 
    c(
        '16p.NSC.DUP',
        '16p.NSC.DEL',
        '16p.NSC.WT',
        '16p.iN.DUP',
        '16p.iN.DEL',
        '16p.iN.WT'
    )
SAMPLE_GROUP_COMPARISONS <- 
    tribble(
        ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        # '16p.iN.DUP',       '16p.iN.DEL', 
        # '16p.NSC.DUP',      '16p.NSC.DEL',
        # '16p.NSC.DUP',      '16p.iN.DUP',
        # '16p.NSC.DEL',      '16p.iN.DEL',
        # '16p.NSC.WT',       '16p.iN.WT',
        # '16p.iN.DUP',       '16p.iN.WT',  
        # '16p.iN.DEL',       '16p.iN.WT',  
        '16p.NSC.DUP',      '16p.NSC.WT',
        '16p.NSC.DEL',      '16p.NSC.WT'
    )
# have consistent colors for genotypes across figures
GENOTYPE_COLORS <- 
    c(
        'WT'='#b3b3ff',
        'DEL'='#ff0000',
        'DUP'='#0000ff'
    )
CELLTYPE_COLORS <- 
    c(
        'NSC'='#b3b3ff',
        'iN'='#ff0000'
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
# Match ideal smoothing param to specified resolution based on this
# https://github.com/TaoYang-dev/hicrep?tab=readme-ov-file#hicrep-parameters
RESOLUTION_IDEAL_H <- 
    tribble(
       ~resolution, ~Ideal_H, 
             10000,       20,
             25000,       10,
             40000,        5,
             50000,        5,
            100000,        3,
            500000,        2,
            500000,        1,
           1000000,        1,
           1000000,        0
    )
# Genomic region boundaries to refernce during analysis plotting
# 16p deleteion and telomere are to replicate Fig 4B in this paper:
# https://www.nature.com/articles/s41588-022-01203-y
GENOMIC_REGIONS <- 
	tribble(
    ~region.group, ~region,                    ~region.chr, ~region.start, ~region.end,
    "16p",         "16p11.2 CNV",              "chr16",          29488679,    30188679,
    "16p",         "chr16p.telomere",          "chr16",                 0,     5149999,
    "16p",         "chr16p11.2",               "chr16",          24288679,    30188679,
    "16p",         "chr16p",                   "chr16",                 0,    36800000,
    "16p",         "chr16q",                   "chr16",          36800001,    90338345,
    "16p",         "chr16",                    "chr16",                 0,    90338345,
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
	) %>%	
    mutate(
        region.UCSC=glue("{region.chr}:{region.start}-{region.end}"),
        region.dist=region.end - region.start
    )
