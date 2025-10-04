###################################################
# Dependencies 
###################################################
library(here)
here::i_am('scripts/run.multiHiCCompare.R')
BASE_DIR <- here()
SCRIPT_DIR <- file.path(BASE_DIR, 'scripts')
source(file.path(SCRIPT_DIR, 'locations.R'))
source(file.path(SCRIPT_DIR, 'constants.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
source(file.path(SCRIPT_DIR, 'utils.multiHiCCompare.R'))
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(purrr)
###################################################
# Load data + Set up comparisons to compute 
###################################################
# GRanges object with Centro/Telomere regions to filter
data('hg38_cyto') 
sample.metadata.df <- 
    load_sample_metadata()
# Set up all combinations of parameters and sample groups to test
covariates.df <- 
    sample.metadata.df %>% 
    select(
        SampleID,
        FlowcellID
    ) %>%
    rename('Batch'=FlowcellID) %>% 
    mutate(Batch=as.factor(Batch))
# All combinations of analysis hyper-params 
hyper.params.df <- 
    expand_grid(
        # zero.p=c(0.6, 0.8),
        zero.p=c(0.8),
        A.min=c(5)
    )
# All pairs of sample groups to compare for differential contacts
comparisons.df <- 
    set_up_sample_comparisons() %>%
# For every group of samples + parameter combination run multiHiCCompare and cache the results.
###################################################
# Run multiHiCCompare on all comparisons
###################################################
# parallelizing params
library(BiocParallel)
numCores <- length(availableWorkers())
message(glue('Using {numCores} workers for parallelism'))
register(MulticoreParam(workers=numCores / 2), default=TRUE)
plan(multisession, workers=numCores / 2)
# Run multiHiCComapre on all comparisons
comparisons.df %>% 
    arrange(desc(chr)) %>% 
    run_all_multiHiCCompare(    
        hyper.params.df=hyper.params.df,
        covariates.df=covariates.df,
        # covariates.df=NULL,
        sample_group_priority_fnc=sample_group_priority_fnc_16p,
        # force_redo=TRUE,
        remove.regions=hg38_cyto,
        p.method='fdr'
    )
# bash one-liner to list results
# find results/multiHiCCompare/results/merged_FALSE/ -type f -name "*-multiHiCCompare.tsv" -exec wc -l {} \; | tr '[/ ]' '\t' |  sed -e "s/-multiHiCCompre.tsv//" | cut -f1,8- |

