###################################################
# Depdendencies
###################################################

###################################################
# Generate Cooltools Results
###################################################
generate_cooltools_calling_cmd <- function(
    threads,
    normalization,
    resolution,
    MatrixID,
    mcool.filepath,
    track.type,
    phasing.track.filepath,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("track.type_{track.type}"),
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}")
        )
    # Create filepaths
    input.filepath     <- glue("{mcool.filepath}::resolutions/{resolution}")
    compartment.prefix <- glue("{output_dir}/{MatrixID}-")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd       <- glue("mkdir -p {output_dir}")
    compartment.cmd <- glue("cooltools eigs-cis --phasing-track {phasing.track.filepath} --n-eigs 3 {weight_flag}  -o {compartment.prefix} {input.filepath}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=glue("{compartment.prefix}.cis.vecs.tsv"),
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    compartment.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_cooltools_calling_cmds <- function(
    hyper.params.df,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # list all hyper-params and corresponding input files (i.e. phasing track files)
    # to call compartments with
    GENOME_TRACK_FILES_DIR %>%
    parse_results_filelist(suffix='-genome.track.tsv') %>%
    dplyr::rename('phasing.track.filepath'=filepath) %>% 
    inner_join(
        hyper.params.df,
        by=join_by(track.type, resolution)
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    cross_join(
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath),
    ) %>% 
    mutate(output_dir=file.path(COMPARTMENTS_RESULTS_DIR, glue("Comp.method_{Comp.method}"))) %>% 
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(Comp.method, ...) {
                        case_when(
                            Comp.method == 'cooltools' ~ generate_cooltools_calling_cmd(...),
                            .unmatched='error'
                        )
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    # Only include cmds generating outputfiles that dont exist
    {
        if (!force_redo) {
            filter(., !file.exists(output.filepath))
        } else {
            .
        }
    }
}

###################################################
# Load cooltools results
###################################################
quantize_track <- function(
    scores,
    nbins=4){
    cut(
        x=scores,
        breaks=
            quantile(
                E1,
                probs=
                    c(
                        seq(0.0, 0.5, length.out=nbins+1),
                        seq(0.5, 1.0, length.out=nbins+1),
                    ) %>%
                    unique()
                na.rm=TRUE
            ),
        labels=
            c(
              paste
            )
    )
}

load_cooltools_compartment_results <- function(
    filepath,
    ...){
    # filepath=tmp$filepath[[1]]
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
        {.} -> e1.df; e1.df
    e1.df %>% 
    # classify bins by first eigenvector
    mutate(
        compartment.n2=
            case_when(
                E1 >  0   ~ 'A',
                E1 <  0   ~ 'B',
                E1 == 0   ~ '0',
                is.na(E1) ~ NA,
                TRUE      ~ as.character(E1)
            )
    ) %>%
    group_by(chrom, compartment.n2) %>% 
    mutate(
        across(
            .fn=,
            .names=''
        )
    )
        compartment.n6=
            cut(
                x=E1,
                breaks=
                    quantile(
                        E1,

                        probs=c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0),
                        na.rm=TRUE
                    ),
                labels=
                    c(
                        'Strong A',
                        'Med A',
                        'Weak A',
                        # 'Near 0',
                        'Weak B',
                        'Med B',
                        'Strong B'
                    )
            )
        compartment.n20=
            cut(
                x=E1,
                breaks=
                    quantile(
                        E1,
                        probs=seq(0, 1, 0.05),
                        na.rm=TRUE
                    ),
                labels=
                    c(
                    )
            )
    ) %>% 
    arrange(chrom, start, end) %>% 
    group_by(chrom) %>% 
    mutate(
        across(
            .cols=starts_with('compartment.'),
            .fn=\(x, ...) glue('{x}->{lead(x, n=1L)}'), # need the ... for dplyr reasons idk
            .names='switch.{str_remove(.col, "^compartment.")}'
        )
    ) %>% 
    ungroup() %>% 
        count(compartment.n6, switch.n6) %>% print(n=Inf)
    select(
        chrom, start, end,
        E1,
        starts_with('compartment.', 'switch.')
    ) %>% 
    filter(!is.na(E1)) %>% 
    dplyr::rename('chr'=chrom)
}

load_all_cooltools_compartment_results <- function(resolutions=NULL){
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix='-.cis.vecs.tsv'
    ) %>%  
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    mutate(
        compartments=
            future_pmap(
                 .l=.,
                 .f=load_cooltools_compartment_results,
                 .progress=TRUE
            )
    ) %>%
    unnest(compartments) %>% 
    select(-c(filepath)) 
}

post_process_cooltools_compartment_results <- function(results.df){
    results.df %>%
    mutate(
        compartment.binary=
            factor(
                compartment.binary,
                levels=
                    c(
                        'A',
                        'B'
                    )
            ),
        compartment.quantile=
            factor(
                compartment.quantile,
                levels=
                    c(
                        'Strong A',
                        '  Weak A',
                        '  Near 0',
                        '  Weak B',
                        'Strong B'
                    )
            ),
        compartment.switch=
            factor(
                compartment.switch,
                levels=
                    c(
                        'A->A'
                        'B->B'
                        'A->B'
                        'B->A'
                    )
            )
    )
}

