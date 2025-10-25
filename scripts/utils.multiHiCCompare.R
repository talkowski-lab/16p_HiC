###################################################
# Depdendencies
###################################################
library(tidyverse)
library(magrittr)
library(glue)
library(multiHiCcompare)
library(BiocParallel)
library(ggplot2)
library(viridis)
library(cowplot)
library(gtable)
library(furrr)

###################################################
# Generate resutls
###################################################
set_up_sample_comparisons <- function(comparison.groups){
    # get info + filepaths for all contact matrices
    load_mcool_files(
        return_metadata_only=TRUE,
        keep_metadata_columns=FALSE
    ) %>% 
    get_min_resolution_per_matrix() %>% 
    distinct() %>% 
    # Now group samples by condition, 
    filter(!isMerged) %>% 
    nest(samples.df=-c(isMerged)) %>% 
    cross_join(comparison.groups) %>%
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.P1.Pattern) ~ Sample.Group.P1,
                        str_detect(SampleID, Sample.Group.P2.Pattern) ~ Sample.Group.P2,
                        TRUE ~ NA
                    )
            ) %>%
            filter(!is.na(Sample.Group)) %>% 
            list()
    ) %>%
    # minimum and max resoltion of all individual matrices per comparison
    mutate(
        resolution.min=min(samples.df$resolution),
        resolution.max=max(samples.df$resolution)
    ) %>% 
    # Now list every comparison at every resolution that is either a min or max for 1 comparison
    ungroup() %>%
    mutate(resolution=list(unique(c(resolution.min, resolution.max)))) %>%
    unnest(resolution) %>% 
    mutate(
        resolution.type=
            case_when(
                resolution == resolution.max ~ 'max',
                resolution == resolution.min ~ 'min',
                TRUE                         ~ NA
            )
    ) %>% 
    select(-c(resolution.min, resolution.max))
}

sample_group_priority_fnc_16p <- function(Sample.Group){
    case_when(
        Sample.Group == '16p.NSC.DUP' ~ 1,  # always numerator in FCs
        Sample.Group == '16p.NSC.DEL' ~ 2,
        Sample.Group == '16p.NSC.WT'  ~ 3,
        Sample.Group == '16p.iN.DUP'  ~ 4,
        Sample.Group == '16p.iN.DEL'  ~ 5,
        Sample.Group == '16p.iN.WT'   ~ 6,
        TRUE                          ~ Inf
    )
}

sample_group_priority_fnc_NIPBLWAPL <- function(Sample.Group){
    # FC is determined by the edger::exactTest() function called in multiHiCCompare
    # https://github.com/dozmorovlab/multiHiCcompare/blob/dcfe4aaa8eaef45e203f3d7f806232bb613d2c9b/R/glm.R#L69
    # According to the docs for exactTest()
    # """Note that the first group listed in the pair is the baseline for the comparison—so if the pair is c("A","B") then the comparison is B - A, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change)."""
    # So with the factor level that comes FIRST is used as the baseline i.e. DENOMINATOR
    # https://www.quantargo.com/help/r/latest/packages/edgeR/NEWS/exactTest

    # so for NIPBL.DEL vs NIPBLWT we want to force NIPBL.DEL to be numerator therefore we make it 
    # the SECOND factor level i.e. have  a larger priority number 
    # i.e. this works when the priority of NIPBL.DEL > NIPBL.WT 
    # tmp.df <- tibble(Sample.Group=c(rep('NIPBL.iN.DEL', 3), rep('NIPBL.iN.WT', 3))) %>% 
    # mutate(group.priority=sample_group_priority_fnc_NIPBLWAPL(Sample.Group), Sample.Group=fct_reorder(Sample.Group, group.priority, .desc=TRUE)) %>%
    # arrange(Sample.Group)
    # tmp.df %>% select(Sample.Group, group.priority)
    # tmp.df %>% pull(Sample.Group)

    case_when(
        grepl(  'All.iN.DEL', Sample.Group) ~ 1, # always numerator since always last factor level
        grepl( 'WAPL.iN.DEL', Sample.Group) ~ 2,
        grepl('NIPBL.iN.DEL', Sample.Group) ~ 3,
        grepl(  'All.iN.WT',  Sample.Group) ~ 4,
        grepl( 'WAPL.iN.WT',  Sample.Group) ~ 5,
        grepl('NIPBL.iN.WT',  Sample.Group) ~ 6,
        TRUE                                ~ Inf
    )
}

handle_covariates <- function(
    samples.df,
    covariates.df,
    effect.col='Sample.Group',
    sampleID.col='SampleID'){
    # Immediate check 
    if (is.null(covariates.df)) {
        covariates <- NULL
        design.matrix <- NULL
        message('No covaraites to use')
    } else {
        # List all covariates for this set of samples
        covariates <- 
            samples.df %>% 
            select(all_of(c(sampleID.col, effect.col))) %>% 
            left_join(
                covariates.df,
                by=sampleID.col
            ) %>%
            select(-all_of(sampleID.col))
        # Check that batch variables are not uniform 
        covariate_level_sizes <- 
            covariates %>% 
            select(-c(all_of(effect.col))) %>% 
            pivot_longer(everything(), names_to='covariate', values_to='value') %>%
            distinct() %>% 
            count(covariate) %>% 
            deframe()
        # If all covariates are uniform, return null, no design matrix
        if (all(covariate_level_sizes == 1)) {
            covariates <- NULL
            design.matrix <- NULL
            message('Covariates are uniform, ignoring')
        } else {
            # List all uniform covariates (uninformative)
            non.uniform.covariates <- 
                covariate_level_sizes[covariate_level_sizes != 1] %>% names()
            # Make formula for GLM
            covariate.names <- 
                covariates %>% 
                select(-c(all_of(effect.col))) %>% 
                colnames() %>%
                {.[. %in% non.uniform.covariates]}
            contrast <- 
                c(effect.col, covariate.names) %>%
                paste(collapse=" + ") %>% 
                sprintf('~ %s', .) %>% 
                formula()
            design.matrix <- model.matrix(contrast, covariates)
        }
    }
    # return final design matrix
    return(
        list(
            covariates=covariates,
            design.matrix=design.matrix
        )
    )
}

run_multiHiCCompare <- function(
    samples.df,
    sample_group_priority_fnc,
    zero.p,
    A.min,
    resolution,
    range1,
    range2,
    remove.regions,
    md_plot_file,
    covariates.df=NULL,
    frac.cutoff=0.8,
    effect.col='Sample.Group',
    p.method='fdr',
    ...){
# row_index=72 * 5 + 16; samples.df=tmp$samples.df[[row_index]]; resolution=tmp$resolution[[row_index]]; range1=tmp$range1[[row_index]]; range2=tmp$range2[[row_index]]; md_plot_file=tmp$md_plot_file[[row_index]]; remove.regions=hg38_cyto; p.method='fdr'; effect.col='Sample.Group'; zero.p=tmp$zero.p[[row_index]]; A.min=tmp$A.min[[row_index]]; frac.cutoff=0.8; samples.df; sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL
    # Handle covariates if specified
    design.info <- 
        handle_covariates(
            samples.df,
            covariates.df,
            effect.col
        )
    # Get sample groups, ensure consistent num/denom for fc estimates
    # function with more priority (smaller number) will be numerator
    # should be no ties, but ties are broken alphabetically
    samples.df <- 
        samples.df %>%
        mutate(
            group.priority=sample_group_priority_fnc(!!sym(effect.col)),
            !!effect.col :=
                fct_reorder(
                    !!sym(effect.col),
                    group.priority,
                    .desc=TRUE
                )
        ) %>%
        arrange(!!effect.col)
    # Load all contacts for samples + regions
    samples.contacts <- 
        samples.df %>%
        select(-c(resolution)) %>% 
        pmap(
            .l=.,
            .f=load_mcool_file,
            resolution=resolution,
            range1=range1,
            range2=range2,
            normalization="NONE",
            .progress=FALSE
        )
    # Get all bin pairs that are detected in > frac.cutoff fraction of samples e.g. 
    # frac.cutoff=0.8 -> only test contacts detected in > 80% of all samples
    common.bin.pairs <- 
        samples.contacts %>%
        bind_rows(.id='index') %>%
        select(chr, range1, range2) %>% 
        count(chr, range1, range2) %>% 
        filter(n >= max(n) * frac.cutoff) %>%
        select(chr, range1, range2)
    message(glue('Testing {nrow(common.bin.pairs)} bin-pairs for DAC'))
    # Now only subset to commonly found contacts 
    samples.contacts <- 
        samples.contacts %>%
        lapply(
            function(df) {
                inner_join(
                    df,
                    common.bin.pairs,
                    by=join_by(chr, range1, range2)
                )
            }
        )
    # Make experiment object with relevant data+parameters
    make_hicexp(
        data_list=samples.contacts,
        groups=samples.df %>% pull(!!effect.col),
        covariates=design.info$covariates,
        zero.p=zero.p,
        A.min=A.min,
        remove.regions=remove.regions,
        remove_zeroes=FALSE,
        filter=TRUE
    ) %>% 
    # Normalize hic data, use cyclic loess with automatically calculated span
    cyclic_loess(
        verbose=TRUE,
        parallel=TRUE,
        span=NA
    ) %T>% 
    # Plot normalized IFs for all sample pairs
    {
        dir.create(
            dirname(md_plot_file),
            showWarnings=FALSE,
            recursive=TRUE
        )
        pdf(
            md_plot_file,
            height=nrow(.@metadata) * 2,
            width=6
        )
        MD_hicexp(
            .,
            pcol=2
        )
        dev.off()
    } %>%
    # Preform differential testing on contacts (bin-pairs)
    # Handle covariate information
    {
        if (is.null(design.info$design.matrix)) {
            hic_exactTest(
                .,
                p.method=p.method,
                parallel=TRUE
            )
        } else {
            hic_glm(
                .,
                design=design.info$design.matrix,
                coef=2,  # Sample.Group (only correct if 2 levels), 1 would be intercept
                method="QLFTest",
                p.method=p.method,
                parallel=TRUE
            )
        }
    } %>% 
    # Get differential results
    results() %>% 
    as_tibble() %>%
    arrange(p.adj)
}

run_all_multiHiCCompare <- function(
    comparisons.df,
    hyper.params.df,
    sample_group_priority_fnc,
    force_redo=FALSE,
    covariates.df=NULL,
    chromosomes=CHROMOSOMES,
    ...){
    # chromosomes=CHROMOSOMES; covariates.df=NULL; force_redo=TRUE; remove.regions=hg38_cyto; sample_group_priority_fnc=sample_group_priority_fnc_NIPBLWAPL; p.method='fdr'
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
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
                Sample.Group.P1.Priority > Sample.Group.P2.Priority ~ Sample.Group.P1,
                Sample.Group.P1.Priority < Sample.Group.P2.Priority ~ Sample.Group.P2,
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
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                MULTIHICCOMPARE_DIR,
                'results',
                glue('merged_{isMerged}'),
                glue('zero.p_{zero.p}'),
                glue('A.min_{A.min}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        md_plot_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-MD.plot.pdf')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-multiHiCCompare.tsv')
            )
    ) %>% 
    arrange(resolution) %>% 
        # {.} -> tmp
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_multiHiCCompare,
                    sample_group_priority_fnc=sample_group_priority_fnc,
                    covariates.df=covariates.df,
                    # all columns also passed as input arguments to run_multiHiCCompare() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this function
        .progress=TRUE
    )
}

###################################################
# Load resutls
###################################################
load_and_correct_multiHiCCompare_results <- function(
    filepaths,
    nom.threshold,
    fdr.threshold,
    gw.fdr.threshold,
    ...){
    filepaths %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    bind_rows() %>%
    mutate(p.adj.gw=p.adjust(p.value, method='BH')) %>% 
    filter(
        p.adj.gw < gw.fdr.threshold,
        p.adj    < fdr.threshold,
        p.value  < nom.threshold
    )
}

load_all_multiHiCCompare_results <- function(
    resolutions=NULL,
    comparisons=NULL,
    gw.fdr.threshold=1,
    fdr.threshold=1,
    nom.threshold=1,
    file_suffix='-multiHiCCompare.tsv',
    ...){
    # comparisons=NULL; resolutions=NULL; gw.fdr.threshold=1; fdr.threshold=1; nom.threshold=0.05; file_suffix='-multiHiCCompare.tsv'
    # Load all results
    parse_results_filelist(
        input_dir=file.path(MULTIHICCOMPARE_DIR, 'results'),
        suffix=file_suffix,
        filename.column.name='pair.name',
        param_delim='_',
    ) %>%
    # Split title into pair of groups ordered by numerator/denominator
    mutate(pair.name=str_remove(pair.name, file_suffix)) %>% 
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=c('Sample.Group.A', 'Sample.Group.B')
    ) %>% 
    rename(
        'Sample.Group.Numerator'=Sample.Group.A,
        'Sample.Group.Denominator'=Sample.Group.B
    ) %>% 
    # Filter relevant results sets
    mutate(
        comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}'),
        resolution=scale_numbers(resolution, force_numeric=TRUE)
    ) %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(comparisons)) {
            filter(., any(grepl(comparisons,  comparison)))
        } else {
            .
        }
    } %>% 
    # Load all results + correct pvalues genome wide per comparison
    group_by(across(-c(filepath, region))) %>% 
    summarize(filepaths=list(c(filepath))) %>% 
    # filter(resolution == 50000) %>% 
    ungroup() %>% 
    mutate(
        results=
            # future_pmap(
            pmap(
                .l=.,
                # correct adjusted pvalues genome-wide
                .f=load_and_correct_multiHiCCompare_results,
                # Only load results with specific params
                gw.fdr.threshold=gw.fdr.threshold,
                fdr.threshold=fdr.threshold,
                nom.threshold=nom.threshold,
                .progress=TRUE
            )
    ) %>% 
        {.} %>% select(comparison, filepaths, results)
    unnest(results) %>% 
    select(-c(filepaths, D)) %>% 
    rename(
        'region1.bp'=region1,
        'region2.bp'=region2
    ) %>% 
    mutate(
        chr=rename_chrs(chr),
        resolution=scale_numbers(resolution),
        distance.bp=region2.bp - region1.bp,
        # calculate log of all pvalues + add columns
        across(
            starts_with('p.'),
            .f=~ -log10(.x),
            .names='log.{.col}'
        )
    ) %>%
    unite(
        'bin.pair.idx',
        chr, region1.bp, region2.bp,
        sep='#',
        remove=FALSE
    )
    # left_join(
    #     load_chr_sizes() %>% rename('chr'=Chr),
    #     by='chr'
    # ) %>% 
    # mutate(
    #     region1.bin=region1.bp / resolution,
    #     region2.bin=region2.bp / resolution,
    #     distance.bin=region2.bin - region1.bin,
    #     region1.pct=100 * (region1.bp / chr.total.bp),
    #     region2.pct=100 * (region2.bp / chr.total.bp),
    #     distance.pct=region2.pct - region1.pct
    # )
}

post_process_multiHiCCompare_results_NIPBLWAPL <- function(results.df){
    results.df %>% 
    mutate(
        chr=factor(chr, levels=CHROMOSOMES),
        resolution=
            resolution %>% 
            scale_numbers(force_numeric=TRUE) %>% 
            scale_numbers(),
        comparison=
            factor(
                comparison,
                levels=
                    c(
                        'NIPBL.iN.DEL vs NIPBL.iN.WT',
                         'WAPL.iN.WT vs NIPBL.iN.WT',
                         'WAPL.iN.DEL vs WAPL.iN.WT',
                        'NIPBL.iN.DEL vs WAPL.iN.WT',
                         'WAPL.iN.DEL vs NIPBL.iN.WT',
                        'NIPBL.iN.DEL vs All.iN.WT',
                         'WAPL.iN.DEL vs All.iN.WT'
                    )
            ),
        comparison.type=
            case_when(
                str_detect(comparison, 'NIPBL.* vs NIPBL.*') ~ 'main',
                str_detect(comparison, 'WAPL.* vs WAPL.*')   ~ 'main',
                str_detect(comparison, '.*WT vs .*WT')       ~ 'main',
                str_detect(comparison, '.*vs All.iN.WT')     ~ 'over.edits',
                str_detect(comparison, 'NIPBL.* vs WAPL.*')  ~ 'across.edits',
                str_detect(comparison, 'WAPL.* vs NIPBL.*')  ~ 'across.edits',
                TRUE                                         ~ '???'
            ) %>%
            factor(levels=c('main', 'across.edits', 'over.edits'))
    )
}

post_process_multiHiCCompare_results_16p <- function(results.df){
    results.df %>% 
    mutate(
        chr=factor(chr, levels=CHROMOSOMES),
        resolution=
            resolution %>% 
            scale_numbers(force_numeric=TRUE) %>% 
            scale_numbers(),
        comparison=
            factor(
                comparison,
                levels=
                    c(
                        '16p.NSC.WT vs 16p.iN.WT',
                        # '16p.NSC.DEL vs 16p.iN.DEL',
                        # '16p.NSC.DUP vs 16p.iN.DUP',
                        '16p.NSC.DUP vs 16p.NSC.DEL',
                        '16p.NSC.DUP vs 16p.NSC.WT',
                        '16p.NSC.DEL vs 16p.NSC.WT',
                        '16p.iN.DUP vs 16p.iN.DEL',
                        '16p.iN.DUP vs 16p.iN.WT',
                        '16p.iN.DEL vs 16p.iN.WT'
                    )
            ),
        comparison.type=
            case_when(
                comparison == '16.NSC.DUP vs 16p.NSC.DEL' ~ 'main.NSC',
                comparison == '16.NSC.DUP vs 16p.NSC.WT'  ~ 'main.NSC',
                comparison == '16.NSC.DEL vs 16p.NSC.WT'  ~ 'main.NSC',
                comparison == '16.NSC.WT vs 16p.iN.WT'    ~ 'wt.vs.wt',
                comparison == '16p.NSC.DEL vs 16p.iN.DEL' ~ 'crosstype',
                        on == '16p.NSC.DUP vs 16p.iN.DUP' ~ 'crosstype',
                comparison == '16.iN.DUP vs 16p.iN.DEL'   ~ 'main.iN',
                comparison == '16.iN.DUP vs 16p.iN.WT'    ~ 'main.iN',
                comparison == '16.iN.DEL vs 16p.iN.WT'    ~ 'main.iN',
                TRUE                                      ~ '???'
            ) %>%
            factor(levels=c('main.NSC', 'main.iN', 'crosstype'))
    )
}

