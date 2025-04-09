library(tidyverse)
library(glue)
BASE_PATH='/data/talkowski/Samples/16p_HiC'
# List all available fastq files
glue('{BASE_PATH}/fastq') %>% 
    list.files(pattern='*16p.*.fastq.gz') %>%
    tibble(filename=.) %>%
    mutate(fastq=glue('{BASE_PATH}/fastq/{filename}')) %>%
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
                'SampleName',
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
    mutate(Sample_Name=glue('{Investigator.ID}_{SampleName}')) %>% 
    left_join(
        glue('{BASE_PATH}/misc/sample.metadata.tsv') %>%
        read_csv() %>%
        select(c(Sample_ID, Sample_Name, Date)),
        by=c('Sample_Name')
    ) %>% 
# extract more info from Sample Name
    mutate(
        CellType=
            case_when(
                str_detect(SampleName, 'NSCHiC$') ~ 'NSC',
                str_detect(SampleName, 'iNHiC')   ~ 'iN',
                TRUE ~ 'Unknown'
            ) %>%
            factor(levels=c('NSC', 'iN')),
        Genotype=
            case_when(
                str_detect(SampleName, '^16pWT')  ~ 'WT',
                str_detect(SampleName, '^16pDEL') ~ 'DEL',
                str_detect(SampleName, '^16pDUP') ~ 'DUP',
                TRUE ~ 'Unknown'
            ) %>%
            factor(levels=c('WT', 'DEL', 'DUP')),
    ) %>%
    write_tsv(SAMPLE_METADATA_TSV)
