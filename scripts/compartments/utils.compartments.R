###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(furrr)
library(HiCDOC)

###################################################
# Generate dcHiC Results
###################################################
load_dcHiC_sample_groups <- function(){
    COMPARTMENTS_PREPROCESSED_DIR %>%
    parse_results_filelist(
        suffix='.matrix',
        filename.column.name='SampleID'
    ) %>%
    get_info_from_SampleIDs() %>% 
    filter(isMerged == 'Individual') %>% 
    dplyr::rename('preprocessed.matrix.filepath'=filepath) %>%
    mutate(
        preprocessed.region.filepath=
            str_replace(
                preprocessed.matrix.filepath,
                '.matrix',
                '.abs.bed'
            )
    ) %>% 
    mutate(
        Sample.Group=as.character(glue('{Edit}.{Celltype}.{Genotype}')),
        Sample.Group.copy=str_replace_all(Sample.Group, fixed('.'), '#'),
        SampleID=str_replace_all(SampleID, fixed('.'), '#')
    ) %>% 
    select(
        Celltype,
        contact.type,
        resolution,
        preprocessed.matrix.filepath,
        preprocessed.region.filepath,
        SampleID,
        Sample.Group,
        Sample.Group.copy
    ) %>% 
    nest(
        input.file.contents=
            c(
                preprocessed.matrix.filepath,
                preprocessed.region.filepath,
                SampleID,
                Sample.Group.copy
            )
    )
}

setup_dcHiC_group_comparisons <- function(
    sample.groups.df,
    comparisons.list,
    cols_to_pair,
    ...){
    sample.groups.df %>% 
    get_all_row_combinations(
        df1=filter(., !str_detect(Sample.Group, '.WT')),
        df2=filter(., str_detect(Sample.Group, '.WT')),
        cols_to_pair=cols_to_pair,
        suffixes=c('.Numerator', '.Denominator'),
        keep_self=FALSE
    ) %>% 
    inner_join(
        comparisons.list,
        by=join_by(Sample.Group.Numerator, Sample.Group.Denominator)
    ) %>% 
    # make input file for dcHiC script
    rowwise() %>% 
    mutate(
        total.input.file.contents=
            list(
                bind_rows(
                    input.file.contents.Numerator,
                    input.file.contents.Denominator
                )
            )
    ) %>% 
    select(-c(input.file.contents.Numerator, input.file.contents.Denominator))
}

make_cmd_list <- function(
    resolution,
    dchic.script.filepath,
    genome_name,
    force_redo,
    input.filepath,
    contact.type,
    threads,
    seed,
    ...){
    base.cmd <- 
        c(
            "Rscript {dchic.script.filepath}",
            "--dirovwt {force_redo}",
            "--seed {seed}",
            "--file {input.filepath}",
            "--output_dir {output_dir}"
        )
    # Generate PCA loadings using cis contacts only
    make.pcs.cmd <- 
        c(
            base.cmd,
            "--pcatype {contact.type}",
            "--cthread {threads}",
            "--pthread {threads}"
        )
        # Pick which PCs to use for compartment assignment
    select.pcs.cmd <- 
        c(
            base.cmd,
            "--pcatype select",
            "--genome {genome_name}",
            "--gfolder {hg38_{resolution}_goldenpathData}"
        )
    # Compute differential statistics between the 2 conditions
    analyze.cmd <- 
        c(
            base.cmd,
            "--pcatype analyze",
            "--diffdir {output.dir}"
        )
    # Annotated subcompartments as well
    subcomp.cmd <- 
        c(
            base.cmd,
            "--pcatype subcomp",
            "--diffdir {output.dir}",
            collapse=" "
        )
    # paste commands together
    list(
        "cd {{output.dir}",
        make.pcs.cmd,
        select.pcs.cmd,
        analyze.cmd,
        subcomp.cmd,
        "cd -"
    ) %>%
    lapply(paste, collapse=" ") %>% 
    paste(collapse="; ") %>% 
    glue()
}

###################################################
# Generate HiCDOC Results
###################################################
run_HiCDOC <- function(
    samples.df,
    resolution,
    chrThreshold,    # Remove chromosomes which are too small to be useful.
    repThreshold,    # Filter sparse replicates to rm uninformative reps with few interactions.
    posThreshold,    # Filter positions (bins) which have too few interactions.
    loessSampleSize,
    ...){
    # row_index=2; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; chrThreshold=tmp$chrThreshold[[row_index]]; repThreshold=tmp$repThreshold[[row_index]]; posThreshold=tmp$posThreshold[[row_index]]; loessSampleSize=tmp$loessSampleSize[[row_index]];
    # data(exampleHiCDOCDataSet)
    row.tmp <- 
    # exampleHiCDOCDataSet %>% 
    HiCDOCDataSetFromCool(
        samples.df$filepath,
        replicates=samples.df$conditions,
        # replicates=samples.df$replicates,
        conditions=samples.df$conditions,
        binSize=resolution
    ) %>% 
    filterSmallChromosomes(threshold=chrThreshold) %>% 
    filterSparseReplicates(threshold=repThreshold) %>% 
    # filterSparseReplicates(threshold=0.99) %>% 
    # filterSparseReplicates(threshold=0.005) %>% 
    # filterWeakPositions(threshold=posThreshold) %>% 
    filterWeakPositions(threshold=0.1) %>% 
    normalizeTechnicalBiases() %>% 
    normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
    detectCompartments()
}

save_HiDOC_results <- function(
    output_dir,
    Sample.Group.Numerator,
    Sample.Group.Denominator,
    ...){
    # parameters(row.tmp)
    compartments(row.tmp)
    differences(row.tmp)
    concordances(row.tmp)
}

run_all_HiCDOC <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    ...){
    # sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL; force_redo=FALSE;
    comparisons.df %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    mutate(
        across(
            matches('Sample.Group.(P1|P2)'),
            ~ sample_group_priority_fnc(.x),
            .names='{.col}.Priority'
        )
    ) %>%
    mutate(
        Sample.Group.Numerator=
            case_when(
                Sample.Group.P1.Priority < Sample.Group.P2.Priority ~ Sample.Group.P1,
                Sample.Group.P1.Priority > Sample.Group.P2.Priority ~ Sample.Group.P2,
                TRUE ~ NA
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == Sample.Group.P1 ~ Sample.Group.P2,
                Sample.Group.Numerator == Sample.Group.P2 ~ Sample.Group.P1,
                TRUE ~ NA
            )
    ) %>% 
    select(-c(matches('Sample.Group.(P1|P2)'))) %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        output_dir=
            file.path(
                COMPARTMENTS_DIR,
                'results',
                'method_HiCDOC',
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('chrThreshold_{chrThreshold}'),
                glue('repThreshold_{repThreshold}'),
                glue('posThreshold_{posThreshold}'),
                glue('loessSampleSize_{loessSampleSize}')
            ),
    ) %>% 
    arrange(resolution) %>% 
        {.} -> tmp
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_HiCDOC,
                    sample_group_priority_fnc=sample_group_priority_fnc,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load Compartment Results
###################################################
