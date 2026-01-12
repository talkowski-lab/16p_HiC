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

make_script_file <- function(
    input_filepath,
    total.input.file.contents,
    output_dir,
    script_filepath,
    dchic.script.filepath,
    force_redo,
    genome_name,
    resolution,
    contact.type,
    threads,
    seed,
    ...){
    # make speicific results directory
    dir.create(
        output_dir,
        recursive=TRUE,
        showWarnings=FALSE
    )
    # write dcHiC specific input file to directory
    write_tsv(
        x=total.input.file.contents,
        file=input_filepath,
        col_names=FALSE
    )
    # make cmds to generate results with dcHiC for this comparison
    base.cmd <- 
        c(
            "Rscript {dchic.script.filepath}",
            "--dirovwt {force_redo}",
            "--seed {seed}",
            "--file {input_filepath}",
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
            "--gfolder {DCHIC_REF_DIR}/hg38_{resolution}_goldenpathData"
        )
    # Compute differential statistics between the 2 conditions
    analyze.cmd <- 
        c(
            base.cmd,
            "--pcatype analyze",
            "--diffdir {output_dir}"
        )
    # Annotated subcompartments as well
    subcomp.cmd <- 
        c(
            base.cmd,
            "--pcatype subcomp",
            "--diffdir {output_dir}",
            collapse=" "
        )
    # paste commands together
    list(
        # "cd {output_dir}",
        make.pcs.cmd,
        select.pcs.cmd,
        analyze.cmd,
        subcomp.cmd
        # "cd -"
    ) %>%
    lapply(paste, collapse=" ") %>% 
    paste(collapse="\n") %>% 
    glue() %>%
    as.character() %>% 
    tibble(cmd=.) %>% 
    write_tsv(
        file=script_filepath,
        col_names=FALSE
    )
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
    # row_index=1; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; chrThreshold=tmp$chrThreshold[[row_index]]; repThreshold=tmp$repThreshold[[row_index]]; posThreshold=tmp$posThreshold[[row_index]]; loessSampleSize=tmp$loessSampleSize[[row_index]];
    tmp1 <- 
        HiCDOCDataSetFromCool(
            paths=samples.df$filepath,
            replicates=samples.df$replicates,
            conditions=samples.df$conditions,
            binSize=resolution
        )
    tmp2 <- 
        tmp1 %>% 
        filterSmallChromosomes(threshold=chrThreshold) %>% 
        # filterSparseReplicates(threshold=repThreshold) %>% 
        filterSparseReplicates(threshold=0.005) %>% 
        # filterWeakPositions(threshold=posThreshold) %>% 
        filterWeakPositions(threshold=0.01) %>% 
        normalizeTechnicalBiases() %>% 
        normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
        detectCompartments()
    names(tmp2)
    str(tmp2)
    compartments(tmp2)
    differences(tmp2)
    concordances(tmp2)
    hicdoc.results.obj <- 
        HiCDOCDataSetFromCool(
            samples.df$filepath,
            replicates=samples.df$replicates,
            conditions=samples.df$conditions,
            binSize=resolution
        ) %>% 
        filterSmallChromosomes(threshold=chrThreshold) %>% 
        filterSparseReplicates(threshold=repThreshold) %>% 
        # filterSparseReplicates(threshold=0.99) %>% 
        # filterSparseReplicates(threshold=0.005) %>% 
        filterWeakPositions(threshold=posThreshold) %>% 
        # filterWeakPositions(threshold=0.1) %>% 
        normalizeTechnicalBiases() %>% 
        normalizeDistanceEffect(loessSampleSize=loessSampleSize) %>% 
        detectCompartments()
    compartments(hicdoc.results.obj)
    differences(hicdoc.results.obj)
    concordances(hicdoc.results.obj)
}

run_all_HiCDOC <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    ...){
    sample_group_priority_fnc=sample_group_priority_fnc_Cohesin; force_redo=FALSE;
    # sample_group_priority_fnc=sample_group_priority_fnc_16p; force_redo=FALSE;
    comparisons.df %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    set_foldchange_direction_as_factor(
        sample_group_priority_fnc=sample_group_priority_fnc,
        priority_col_pattern=starts_with('Sample.Group.'),
        label_cols=c('Sample.Group.Left', 'Sample.Group.Right')
    ) %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        output_dir=
            file.path(
                COMPARTMENTS_RESULTS_DIR,
                'method_HiCDOC',
                glue('resolution_{scale_numbers(resolution)}'),
                glue('chrThreshold_{chrThreshold}'),
                glue('repThreshold_{repThreshold}'),
                glue('posThreshold_{posThreshold}'),
                glue('loessSampleSize_{loessSampleSize}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-HiCDOC.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !(file.exists(results_file)))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            dplyr::count(
                .,
                resolution,
                chrThreshold,
                repThreshold,
                posThreshold,
                loessSampleSize,
                Sample.Group.Numerator,
                Sample.Group.Denominator
            ) %>%
            rename_with(~ str_remove(., 'Sample.Group.'))
        )
    } %>%
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
                    ...  # passed from the call run_HiCDOC()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load Compartment Results
###################################################
