library(stringi)
library(furrr)
library(idr2d)
# library(plyranges)
###################################################
# Handle cCRE data
###################################################
load_CTCF_sites <- function(force.redo=FALSE){
    # https://dozmorovlab.github.io/CTCF/
    # annotation object
    check_cached_results(
        results_file=CTCF_SITE_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                AnnotationHub() %>% 
                # subset specific set of CTCF sites predicted in humans by JASPAR 2022
                subset(
                    species == 'Homo sapiens' &
                    genome == 'hg38' &
                    dataprovider == 'JASPAR 2022' &
                    preparerclass == "CTCF"
                ) %>%
                {.[['AH104727']]} %>%
                # make tidy tibble
                as_tibble() %>% 
                dplyr::rename(
                    'length'=width,
                    'chr'=seqnames
                ) %>% 
                select(
                    -c(
                        # db,
                        # db.set,
                        # bioset,
                        # genomic.feature,
                        sequence
                    )
                )
            }
    )
}

load_encode_ccres <- function(force.redo=FALSE){
    check_cached_results(
        results_file=ENCODE_CCRE_ANNOTATIONS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                # Graphical + explicit definitions of of cCRE types
                # https://screen-v4.wenglab.org/about
                # cCRE.descriptions <- 
                #     tribble(
                #         ~cCRE.type,   ~cCRE.description,
                #         'CA',         'Chromatin Accessible Only',
                #         'CA-CTCF',    'Chromatin Accessible with CTCF',
                #         'CA-H3K4me3', 'Chromatin Accessible with H4K4me3',
                #         'CA-TF',      'Chromatin Accessible with TF',
                #         'PLS',        'Promoter-like',
                #         'TF',         'TF Only',
                #         'dELS',       'Distal enhancer-like'
                #         'pELS',       'Proximal enhancer-like'
                #     )
                ccres.df <- 
                    file.path(ENCODE_CCRE_DIR, 'GRCh38-cCREs.bed') %>%
                    read_tsv(
                        col_names=
                            c(
                                'chr',
                                'start', 'end',
                                # 'cCREID.P1', 'cCREID.P2',
                                'cCREID', 'cCREID.2',
                                'cCRE.type'
                            )
                    ) %>% 
                    select(-c(cCREID.2))
                    # left_join(cCRE.descriptions, by='cCRE.type')
            }
    )
}

split_cCRE_annotations <- function(cCRES.df){
    # Split large file into 1 filer per cCRE type 
    # easier to count per type with bedtools this way (multiple databases)
    cCRES.df %>% 
    nest(positions=-c(cCRE.type)) %>% 
    mutate(
        output_dir=
            file.path(
                CCRES_BED_DIR,
                glue('cCRE.type_{cCRE.type}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{cCRE.type}_ENCODE.cCREs.bed')
            )
    ) %>%
    pwalk(
        .l=.,
        .f=
            function(positions.df, results_file, ...){
                dir.create(dirname(results_file), showWarnings=FALSE, recursive=TRUE)
                write_tsv(
                    positions.df,
                    results_file,
                    col_names=FALSE
                )
            }
    )
}

list_all_cCRES_bed_files <- function(){
    CCRES_BED_DIR %>%
    parse_results_filelist(
        suffix='.bed',
        filename.column.name='tmp',
        parse_filepath_to_columns=TRUE
    ) %>%
    dplyr::rename('cCRE.filepath'=filepath) %>% 
    select(-c(tmp))
}

###################################################
# Make BED Files
###################################################
# BED utils
make_region_bed_files <- function(regions.df){
}

make_anchors_bed_files <- function(regions.df){
}

write_all_bed_files <- function(
    features.df,
    force_redo=FALSE){
    features.df %>% 
    {
        if (!force_redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %>% 
    # pmap(
    future_pwalk(
        .l=.,
        .f=
            function(features, results_file, ...){
                dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
                features %>%
                arrange(across(colnames(.)[1:3])) %>% 
                write_tsv(
                    str_replace_all(results_file, ' ' , '.'),
                    col_names=TRUE
                )
            },
        .progress=TRUE
    )
}
# TADs
format_TADs_for_bed_files <- function(TADs.df){
    TADs.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    nest(
        features=
            c(
                `# chr`,
                start, 
                end,
                starts_with('TAD.'),
                -TAD.params
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                TAD_BED_FILES_DIR,
                param_dir,
                'region_TAD.spans',
                glue('{Sample.Group}-TADs.bed')
            )
    )
}

format_TAD_anchors_for_bed_files <- function(TADs.df){
    TADs.df %>% 
    # The end is actually the end of the last bin, transformt so its the start of the last bin inside the TAD
    mutate(end=end - resolution) %>% 
    pivot_TADs_to_boundaries() %>% 
    select(-c(TAD.index)) %>% 
    mutate(boundary.end=boundary + resolution) %>% 
    dplyr::rename(
        "boundary.start"=boundary,
        "boundary.score"=score,
        '# chr'=chr
    ) %>% 
    nest(
        features=
            c(
                `# chr`,
                boundary.start, 
                boundary.end,
                boundary.side,
                boundary.score,
                TAD.bins,
                TAD.length

            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                TAD_BED_FILES_DIR,
                param_dir,
                'region_TAD.boundaries',
                glue('{Sample.Group}-TAD.boundaries.bed')
            )
    )
}
# Loops
make_loop_span_bed_files <- function(loops.df) {
    loops.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    mutate(anchor.right=anchor.right + resolution) %>%
    nest(
        features=
            c(
                `# chr`,
                anchor.left,
                anchor.right,
                length,
                count,
                enrichment,
                log10.qval
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'region_loop.spans',
                glue('{SampleID}-loops.bed')
            )
    )
}

make_loop_bedpe_files <- function(loops) {
    loops.df %>% 
    dplyr::rename(
        'anchor.left.start'=anchor.left,
        'anchor.right.start'=anchor.right,
        '# chr.left'=chr
    ) %>% 
    mutate(
        chr.right=`# chr.left`,
        anchor.left.end=anchor.left.start   + resolution,
        anchor.right.end=anchor.right.start + resolution
    ) %>% 
    nest(
        features=
            c(
                `# chr.left`, anchor.left.start,  anchor.left.end,
                   chr.right, anchor.right.start, anchor.right.end,
                length,
                count,
                enrichment,
                log10.qval
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'region_loops',
                glue('{SampleID}-loops.bedpe')
            )
    )
}

make_loop_anchor_bed_files <- function(loop.valency.df) {
    loop.valency.df %>% 
    mutate(
        param_dir=
            file.path(
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('region_{chr}')
            )
    ) %>% 
    # pivot stats so every row represents a unique loop
    pivot_wider(
        names_from=feature,
        names_glue="loops_{feature}_{.value}",
        values_from=c(min, mean, median, max)
    ) %>% 
    # Prepare for writing to bed files
    dplyr::rename(
        'anchor.start'=anchor.position,
        '# chr'=chr
    ) %>% 
    mutate(anchor.end=anchor.start + resolution) %>% 
    # add so bed file can incluide column titles in the header
    nest(
        features=
            c(
                `# chr`,
                anchor.start,
                anchor.end,
                valency,
                starts_with('loops_')
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'region_loop.anchors',
                glue('{SampleID}-loop.anchors.bed')
            )
    )
}

make_loop_nesting_bed_files <- function(loop.nesting.df) {
    loop.nesting.df %>% 
    mutate(
        param_dir=
            file.path(
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('region_{chr}')
            )
    ) %>% 
    # Prepare for writing to bed files
    dplyr::rename('# chr'=chr) %>% 
    mutate(nest.end=nest.end + resolution) %>% 
    # add so bed file can incluide column titles in the header
    nest(
        features=
            c(
                `# chr`,
                nest.start,
                nest.end,
                nesting.lvl,
                starts_with('metric.')
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'region_loop.nesting',
                glue('{SampleID}-loop.nesting.bed')
            )
    )
}
# DIRs
make_DIR_span_bed_files <- function(DIRs.df){
    DIRs.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    nest(
        features=
            c(
                `# chr`,
                region1,
                region2,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'region_DIR.spans',
                glue('{Sample.Group}-DIRs.bed')
            )
    )
}

make_DIR_bedpe_files <- function(DIRs.df) {
    DIRs.df %>% 
    dplyr::rename(
        'region1.start'=region1,
        'region2.start'=region2,
        '# chr.left'=chr
    ) %>% 
    mutate(
        chr.right=`# chr.left`,
        region1.end=region1.start + resolution,
        region2.end=region2.start + resolution
    ) %>% 
    nest(
        features=
            c(
                `# chr.left`, region1.start, region1.end,
                   chr.right, region2.start, region2.end,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'region_DIRs',
                glue('{Sample.Group}-DIRs.bedpe')
            )
    )
}

make_DIR_anchor_bed_files <- function(DIRs.df) {
    DIRs.df %>% 
    dplyr::rename(
        '# chr'=chr,
        'anchor.left'=region1,
        'anchor.right'=region2,
    ) %>% 
    pivot_longer(
        c(anchor.left, anchor.right),
        names_to='anchor.side',
        values_to='anchor.start'
    ) %>%
    mutate(anchor.end=anchor.start + resolution) %>% 
    nest(
        features=
            c(
                `# chr`,
                anchor.start,
                anchor.end,
                anchor.side,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'region_DIR.anchors',
                glue('{Sample.Group}-DIR.anchors.bed')
            )
    )
}

###################################################
# Build BED Commands
###################################################
list_bed_files <- function(
    bed_dir,
    suffix='-.*.bed',
    filename.column.name='SampleID',
    ...){
    bed_dir %>% 
    parse_results_filelist(
        suffix=suffix,
        filename.column.name=filename.column.name
    ) %>%
    mutate(
        param_dir=dirname(str_remove(filepath, glue('{bed_dir}/'))),
        file.type=str_extract(filepath, '\\.([^\\.]+)$', group=1),
        results_file=
            file.path(
                FUNCTIONAL_ENRICHMENTS_DIR,
                basename(bed_dir),
                param_dir,
                glue('{SampleID}-intersect.cCRE.tsv')
            )
    ) %>%
    select(
        resolution,
        region,
        file.type,
        SampleID,
        param_dir,
        filepath,
        results_file
    )
}

make_bedtools_cmds <- function(
    feature.bed.files.df,
    bed_cmds_file,
    force_redo=FALSE){
    dir.create(dirname(bed_cmds_file), showWarnings=FALSE, recursive=TRUE)
    feature.bed.files.df %>% 
    {
        if (!force_redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %>% 
    # cCRE bed files to count overlaps with for each features
    cross_join(
        list_all_cCRES_bed_files() %>% 
        mutate(cCRE.filepath=str_replace_all(cCRE.filepath, ' ' , '.')) %>% 
        summarize(
            cCRE.filepaths=paste(cCRE.filepath, collapse=' '),
            cCRE.names=paste(cCRE.type, collapse=' '),
        )
    ) %>% 
    # generate bedtools intersect command
    # bedtools intersect 
    #     -a loops.bed    # all loops called on the same chr at the same resolution
    #     -b cCRE_files   # list of all cCRE database files to intersect against
    #     -wao            # keep all overalpping region info
    #     -C              # count how many cCREs of each type overlap each region
    #     -names          # include cCREs names in counting
    mutate(
        # remove spaces in filenames
        across(
            c(filepath, results_file),
            ~ str_replace_all(.x, ' ', '.')
        ),
        mkdir.cmd=
            glue('mkdir -p {dirname(results_file)}'),
        add.column.names.cmd=
            glue("head -1 {filepath} | sed -e 's/^# //' | sed -e 's/$/\\tcCRE.type\\tcCRE.overlaps/' >| {results_file}"),
        bedtools.intersect.cmd=
            glue('bedtools intersect -a {filepath} -b {cCRE.filepaths} -names {cCRE.names} -wao -C >> {results_file}'),
        bed.cmd=
            glue('{mkdir.cmd}; {add.column.names.cmd}; {bedtools.intersect.cmd}')
    ) %>% 
    select(bed.cmd) %>%
    write_tsv(
        # file.path(LOOPS_DIR, 'all.loop.nesting.bedtools.cmds.txt'),
        bed_cmds_file,
        col_names=FALSE
    )
}

###################################################
# Load cCRE enrichment data
###################################################
load_cCRE_overlap_results <- function(filepath){
}

list_all_cCRE_overlap_results <- function(hic.features.df){
    hic.features.df %>% 
    parse_results_filelist(suffix='-intersect.cCRE.tsv')
}

load_all_cCRE_overlap_results <- function(hic.features.df){
    hic.features.df %>% 
    list_all_cCRE_overlap_results() %>% 
    mutate(
        cCRE.overlaps=
            # pmap(
            future_pmap(
                .,
                load_cCRE_overlap_results,
                .progress=TRUE
            )
    ) %>%
        # {.} -> tmp; tmp
        # tmp %>% select(SampleID, loops)
    unnest(cCRE.overlaps) %>% 
    select(
        -c(
            filepath
        )
    )
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

