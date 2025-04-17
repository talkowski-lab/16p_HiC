library(tidyverse)
library(glue)
library(here)
source(file.path(here(), 'scripts/locations.R'), chdir=TRUE)
print(BASE_DIR) # /data/talkowski/Samples/16p_HiC
SAMPLE_INFO_REGEX <- 
    sprintf(
        "(?<Edit>%s)(?<Genotype>%s)(?<SampleNumber>[a-zA-Z0-9]+)(?<Celltype>%s)HiC",
        paste(EDITS, collapse="|"),
        paste(GENOTYPES, collapse="|"),
        paste(CELLTYPES, collapse="|")
    )
# List all available fastq files
FASTQ_DIR %>% 
list.files(pattern='*16p.*.fastq.gz') %>%
tibble(filename=.) %>%
mutate(fastq=glue('{BASE_DIR}/fastq/{filename}')) %>%
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
            'Sample.Info',
            NA,
            NA,
            'Fastq.Side',
            NA
        )
) %>%
# 1 row per sample, with R1,R2 files as columns
pivot_wider(
    names_from=Fastq.Side,
    names_prefix='Fastq.',
    values_from=fastq 
) %>% 
# Get date experiments were submitted
mutate(Sample_Name=glue('{Investigator.ID}_{Sample.Info}')) %>%
left_join(
    glue('{BASE_DIR}/misc/SamplesSheet.tsv') %>%
    read_csv(show_col_types=FALSE) %>%
    select(c(Sample_Name, Date)),
    by=join_by(Sample_Name)
) %>% 
# extract data for each sample
mutate(
    Parsed.Sample.Info=
        Sample.Info %>%
        str_match(SAMPLE_INFO_REGEX) %>% 
        {.[,2:ncol(.)]} %>% 
        as_tibble()
) %>% 
unnest(Parsed.Sample.Info) %>% 
rename(Sample.Number=SampleNumber) %>% 
mutate(
    Edit=factor(Edit, levels=EDITS),
    Genotype=factor(Genotype, GENOTYPES),
    Celltype=factor(Celltype, CELLTYPES)
) %>%
# Convienience/Cleaning
mutate(
    Sample.ID=glue('{Edit}.{Genotype}.{Sample.Number}.{Celltype}'),
    Sample.Name=glue('{Genotype}.{Sample.Number}')
) %>% 
select(-c(Sample_Name)) %>% 
write_tsv(SAMPLE_METADATA_TSV)
print(glue("Sample metadata available at {SAMPLE_METADATA_TSV}"))
