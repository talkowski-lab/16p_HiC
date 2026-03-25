library(stringi)
library(furrr)
library(idr2d)
# library(plyranges)

###################################################
# BED utils
###################################################
make_region_bed_files <- function(regions.df){
}

make_anchors_bed_files <- function(regions.df){
}

list_all_bed_files <- function(suffix){
    BED_FILES_DIR %>%
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix=suffix
    ) %>%
    mutate()
}

make_bedtools_cmds <- function(){
    list_all_bed_files() %>% 
    mutate(
        output_dir=
            file.path(
                ALL_LOOP_NESTING_RESULTS_DIR,
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID}-cCRE.overlaps.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %>% 
    # bedtools intersect 
    #     -a loops.bed    # all loops called on the same chr at the same resolution
    #     -b cCRE_DIR     # dir where each file is all locations of a single type of cCRE
    #     -wao            # keep all overalpping region info
    #     -f 1.0          # only keep cCRES fully within loops
    #     -C              # count how many cCREs of each type overlap each region
    #     -filenames      # include cCREs 
    mutate(bed.cmd=glue('mkdir -p {output_dir} && bedtools intersect -a {filepath} -b {CCRES_BED_DIR} -wao -C -f 1.0 >| {results_file}')) %>% 
    select(bed.cmd) %>%
    write_tsv(
        file.path(LOOPS_DIR, 'all.loop.nesting.bedtools.cmds.txt'),
        col_names=FALSE
    )
}

###################################################
# Functional Enrichment: cCRE 
###################################################
# Functional Enrichment: cCRE ~ TADs
count_cCREs_per_TAD <- function(
    cCREs.df,
    TADs.df, 
    regions.df=NULL){
}

# Functional Enrichment: cCRE ~ Loop 
count_cCREs_per_loop <- function(
    cCREs.df=NULL,
    idr2d.df=NULL){
    # cis-regulatory element annotations (cCREs)
    # loop.totals.df
    # ccre.totals.df
    idr2d.df %>%
    # nest so one set of loop annotations per row
    nest(
        loops=
            c(
                # is.loop.shared,
                diff.value,
                diff.rank,
                IDR,
                anchor.left,
                anchor.right
            )
    ) %>% 
    # match loop sets and cCRE sets by region
    left_join(
        # nest so one set of cCREs annotations per row
        cCREs.df %>%
        nest(
            cCREs=
                c(
                    # cCRE.type,
                    start,
                    end,
                    cCREID
                )
        ),
        by=
            join_by(
                chr,
                region
            )
    ) %>% 
    # For each region, compute all overlaps between any loop and any cCRE
    rowwise() %>% 
    mutate(
        overlaps=
            inner_join(
                loops,
                cCREs,
                by=
                    join_by(
                        within(
                            y$start, y$end,
                            x$anchor.left, x$anchor.right
                        )
                    )
            ) %>%
            list()
    ) %>% 
    ungroup() %>% 
        {.} -> tmp; tmp
        tmp %>% 
            rowwise() %>% mutate(across(c(loops, cCREs, overlaps), ~ nrow(.x))) %>% ungroup() %>%
            dplyr::rename('max.gap'=max.gap.bins.int) %>% 
            select(
                   resolution, max.gap,
                   comparison,
                   region,
                   is.loop.shared, loops,
                   cCRE.type, cCREs,
                   overlaps
            )
        tmp %>% 
            rowwise() %>% mutate(across(c(loops, cCREs, overlaps), ~ nrow(.x))) %>% 
            ungroup() %>% summarize(across(c(loops, cCREs, overlaps), sum))
}

# Functional Enrichment: cCRE ~ Compartment
count_cCREs_per_compartment <- function(
    cCREs.df=NULL,
    TADs.df=NULL, 
    regions.df=NULL){
}
###################################################
# Expression Integration Analysis
###################################################
calc_expr_loop_ztest <- function(
    idr2d.results.df,
    expression.df,
    ...){
    idr2d.results.df %>% 
    inner_join(
        expression.df,
        by=
            join_by(
                SampleID.P1 == SampleID,
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    ) %>% 
    inner_join(
        expression.df,
        suffix=c('.P1', '.P2'),
        by=
            join_by(
                SampleID.P2 == SampleID,
                chr,
                start, end,
                symbol, EnsemblID
            )
    ) %>%
    # for each gene compute pvalue if mean expression is different between conditions
    add_column(
        n.rna.replicates.P1=6,
        n.rna.replicates.P2=6
    ) %>% 
    mutate(
        TPM.se.P1=TPM.sd.P1**2 / n.rna.replicates.P1,
        TPM.se.P2=TPM.sd.P2**2 / n.rna.replicates.P2,
        expr.Z=(TPM.mean.P2 - TPM.mean.P1) / sqrt(TPM.se.P1 + TPM.se.P2),
        expr.p=2 * (1 - pnorm(abs(expr.Z)))
    ) %>%
    # adjust p-values genome-wide
    group_by(
        weight, resolution, kernel,
        resolve.method, metric,
        SampleID.P1, SampleID.P2
    ) %>% 
    mutate(expr.p.adjust=p.adjust(expr.p, method='BH')) %>% 
    ungroup() %>% 
    select(
        -c(
            n.rna.replicates.P1, n.rna.replicates.P2,
            TPM.se.P1, TPM.se.P2,
            # TPM.sd.P1, TPM.sd.P2,
            # TPM.mean.P1, TPM.mean.P2,
            expr.Z, expr.p
        )
    )
}

