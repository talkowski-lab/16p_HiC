library(TADCompare)

###################################################
# Generate TADCompare results
###################################################
TADCompare_load_matrix <- function(
    filepath,
    ...){
    load_mcool_file(
        filepath,
        type='df',
        cis=TRUE,
        ...
    ) %>% 
    select(c(range1, range2, IF))
}

run_TADCompare <- function(
    filepath.Numerator,
    SampleID.Numerator,
    pre_tads.Numerator,
    filepath.Denominator,
    SampleID.Denominator,
    pre_tads.Denominator,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    ...){
    # paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', collapse='; ')
    # chr1 @ 10Kb -> 24896x24896 matrix -> 40Gb is enough
    # Run TADCompare on the 2 matrices being compared
    matrix.numerator <-
        TADCompare_load_matrix(
            filepath.Numerator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    matrix.denominator <-
        TADCompare_load_matrix(
            filepath.Denominator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    pre_tads <- 
        if (is.null(pre_tads.Numerator)) {
            NULL
        } else {
            list(pre_tads.Numerator, pre_tads.Denominator)
        }
    tad.compare.results <- 
        TADCompare(
            matrix.numerator,
            matrix.denominator,
            resolution=resolution,
            z_thresh=z_thresh,
            window_size=window_size,
            gap_thresh=gap_thresh,
            pre_tads=pre_tads
        )
    # Format results to include boundary+gap scores for all bins + differential annotations
    # tad.compare.results$Boundary_Scores %>% as_tibble()
    # tad.compare.results$TAD_Frame %>% as_tibble()
    tad.compare.results$Boundary_Scores %>% 
    as_tibble() %>%
    full_join(
        tad.compare.results$TAD_Frame %>%
        as_tibble() %>% 
        add_column(isTADBoundary=TRUE),
        suffix=c('.All', '.TADs'),
        by=join_by(Boundary)
    ) %>% 
        # {.} -> tcr; tcr
        # tcr %>% count(isTADBoundary, Differential.All, Differential.TADs, Type.All, Type.TADs)
        # tcr %>% 
    mutate(
        isTADBoundary=ifelse(is.na(isTADBoundary), FALSE, isTADBoundary),
        Differential=
            case_when(
                is.na(Differential.TADs) ~ Differential.All,
                TRUE                     ~ Differential.TADs
            ),
        is.Differential=!grepl('Non-Differential', Differential),
        Type=
            case_when(
                is.na(Type.TADs) ~ Type.All,
                TRUE             ~ Type.TADs
            ),
        Enriched.Condition=
            case_when(
                Enriched_In.TADs == 'Matrix 1' ~ SampleID.Numerator,
                Enriched_In.TADs == 'Matrix 2' ~ SampleID.Denominator,
                Enriched_In.All  == 'Matrix 1' ~ SampleID.Numerator,
                Enriched_In.All  == 'Matrix 2' ~ SampleID.Denominator,
                TRUE                      ~ NA
            ),
        TAD_Score1=
            case_when(
                is.na(TAD_Score1.TADs) ~ TAD_Score1.All,
                TRUE                   ~ TAD_Score1.TADs
            ),
        TAD_Score2=
            case_when(
                is.na(TAD_Score2.TADs) ~ TAD_Score2.All,
                TRUE                   ~ TAD_Score2.TADs
            ),
        Gap_Score=
            case_when(
                is.na(Gap_Score.TADs) ~ Gap_Score.All,
                TRUE                  ~ Gap_Score.TADs
            )
    ) %>%
    dplyr::rename(
        'TAD.Score.Numerator'=TAD_Score1,
        'TAD.Score.Denominator'=TAD_Score2,
        'TAD.isDifferential'=Differential,
        'TAD.Difference.Type'=Type
    ) %>% 
    rename_with(~ str_replace_all(.x, '_', '.')) %>% 
    dplyr::select(-c(ends_with('.All'), ends_with('.TADs')))
}

run_all_TADCompare <- function(
    comparisons.df,
    hyper.params.df,
    force_redo=FALSE,
    ...){
    # force_redo=TRUE;
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TADCOMPARE_DIR,
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('TAD.method_{TAD.method}'),
                glue('TAD.params_{TAD.params}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                # glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.Numerator}_vs_{SampleID.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    # filter(chr != 'chrY') %>% 
    arrange(desc(resolution), desc(chr)) %>% 
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
                # z_thresh,
                # window_size,
                # gap_thresh,
                # TAD.params,
                TAD.method, 
                resolution,
                SampleID.Numerator,
                SampleID.Denominator
            )
        )
    } %>%
        # {.} -> tmp; tmp
        # tmp %>% 
    # pmap(
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_TADCompare,
                    # all args also passed as input arguments to run_all*() by pmap
                    ...  # passed from the call to this wrapper()
                )
            },
        ...,  # passed from the call to from run_all_TADCompare
        .progress=TRUE
    )
}

load_TADCompare_results <- function(
    filepath,
    boundaries_only=TRUE,
    ...){
    # row_index=1; filepath=tmp$filepath[[row_index]];
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    {
        if (boundaries_only) {
            filter(., isTADBoundary)
        } else {
            .
        }
    } %>% 
    filter(TAD.Difference.Type != 'Non-Differential')
}

load_and_correct_TADCompare_results <- function(
    filepaths,
    nom.threshold,
    fdr.threshold,
    gw.fdr.threshold,
    ...){
    # filepaths=tmp$filepaths[[1]]
    filepaths %>% 
    mutate(
        results=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                id='tmpID',
                show_col_types=FALSE,
                progress=FALSE
            )
    ) %>% 
    unnest(results) %>% 
        # {.} -> lactr.tmp; lactr.tmp
        # lactr.tmp %>% count(tmpID)
        # lactr.tmp %>% 
    # calculate pvalue from Gap Score (a z-score) calcualted by TADCompare
    # See Section 2.7 here
    # https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00158/full 
    mutate(p.value=2 * (1 - pnorm(abs(Gap.Score)))) %>% 
    # genome-wide fdr adj
    mutate(p.adj.gw=p.adjust(p.value, method='BH')) %>% 
    # chr-wise fdr adjustment
    group_by(tmpID) %>% 
    mutate(p.adj=p.adjust(p.value, method='BH')) %>% 
    ungroup() %>% 
    # filter insiginicant differences
    filter(
        p.adj.gw < gw.fdr.threshold,
        p.adj    < fdr.threshold,
        p.value  < nom.threshold
    ) %>%
    select(-c(tmpID))
}

list_all_TADCompare_results <- function(){
    # Get a list of all results files
    TADCOMPARE_DIR %>% 
    parse_results_filelist(
        suffix='-TADCompare.tsv',
        filename.column.name='pair.name'
    ) %>% 
    # Split title into pair of groups ordered by numerator/denominator
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=c('Sample.Group.Numerator', 'Sample.Group.Denominator')
    )
}

load_all_TADCompare_results <- function(
    nom.threshold,
    fdr.threshold,
    gw.fdr.threshold,
    ...){
    # gw.fdr.threshold=1; fdr.threshold=0.1; nom.threshold=0.05
    list_all_TADCompare_results() %>% 
    # filter(TAD.method != 'cooltools') %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
        # {.} -> tmp
    mutate(
        results=
            # pmap(
            future_pmap(
                .l=.,
                # load_TADCompare_results,
                .f=load_and_correct_TADCompare_results,
                nom.threshold=nom.threshold,
                fdr.threshold=fdr.threshold,
                gw.fdr.threshold=gw.fdr.threshold,
                boundaries_only=TRUE,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

load_correct_count_TADCompare_results <- function(
    filepaths,
    sig.colname='p.adj.gw',
    ...){
    filepaths %>% 
    load_and_correct_TADCompare_results(
        nom.threshold=1,
        fdr.threshold=1,
        gw.fdr.threshold=1
    ) %>% 
    # for each thresh, make binary col if TAD difference meets threshold
    mutate(
        "sig.lvl.{sig.colname} < 1e-15" := .data[[sig.colname]] <  1e-15,
        "sig.lvl.{sig.colname} < 1e-10" := .data[[sig.colname]] <  1e-10,
        "sig.lvl.{sig.colname} < 1e-05" := .data[[sig.colname]] <  1e-5,
        "sig.lvl.{sig.colname} < 0.001" := .data[[sig.colname]] <  1e-3,
        "sig.lvl.{sig.colname} < 0.05 " := .data[[sig.colname]] <  0.05,
        "sig.lvl.{sig.colname} < 0.1  " := .data[[sig.colname]] <  0.10,
        "sig.lvl.N.S."                  := .data[[sig.colname]] >= 0.10
    ) %>% 
    pivot_longer(
        starts_with('sig.lvl.'),
        names_to='sig.lvl',
        names_prefix='sig.lvl.',
        values_to='meet.sig.lvl'
    ) %>% 
    # Inclusively count how many TAD differences meet each thrshold across categories
    # This produces inclusive counts  for each significance threshold i.e. 
    # the number of TAD differences < 0.1 also includes all differences <= 0.01
    filter(meet.sig.lvl) %>% 
    count(
        isTADBoundary,
        TAD.Difference.Type,
        Enriched.Condition,
        region,
        sig.lvl
    ) %>%
    # order significance categories by increasing significance
    mutate(
         sig.lvl=
             factor(
                 sig.lvl,
                 levels=
                     c(
                         'N.S.',
                         as.character(glue('{sig.colname} < 0.1  ')),
                         as.character(glue('{sig.colname} < 0.05 ')),
                         as.character(glue('{sig.colname} < 0.001')),
                         as.character(glue('{sig.colname} < 1e-05')),
                         as.character(glue('{sig.colname} < 1e-10')),
                         as.character(glue('{sig.colname} < 1e-15'))
                     )
             )
     )
}

load_correct_count_all_TADCompare_results <- function(){
    list_all_TADCompare_results() %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
    mutate(
        results=
            future_pmap(
                .l=.,
                .f=load_correct_count_TADCompare_results,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        # 'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

post_process_TADCompare_results <- function(results.df){
    results.df %>%
    # filter(TAD.method != 'cooltools') %>% 
    mutate(
        across(
            c(
                Sample.Group.Numerator,
                Sample.Group.Denominator,
                Enriched.Condition
            ),
            ~ str_remove(.x, '.Merged.Merged')
        ),
        comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}'),
        isBoundary=ifelse(isBoundary, 'TAD', 'Not TAD')
    ) %>% 
    relocate(
        c(
            resolution,
            TAD.params,
            TAD.method,
            Sample.Group.Numerator, Sample.Group.Denominator, 
            comparison, Enriched.Condition,
            chr,
            isBoundary, DifferenceType
        )
    ) %>% 
    select(
        -c(
            z.thresh,
            window.size,
            gap.thresh
        )
    )
}

