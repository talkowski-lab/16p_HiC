###################################################
# Dependencies
###################################################
library(tidyverse)
library(magrittr)
library(tictoc)
library(glue)
library(optparse)
library(future)
library(HiCExperiment)
library(hictkR)

###################################################
# Pairsing and Caching
###################################################
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        message("No results file, will just return data")
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        # Set read/write functions based on filetype
        output_filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (output_filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output_filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {extension}'))
        }
        if (file.exists(results_file) & !force_redo) {
            message(glue('{results_file} exists, not recomputing results'))
            if (return_data) {
                message('Loading cached results...')
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                message(glue('{results_file} exists, recomputing results anyways'))
            } else {
                message(glue('No cached results, generating: {results_file}'))
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    # Or dont save the results, just return the results
    toc()
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

parse_results_filelist <- function(
    input_dir,
    suffix,
    filename.column.name='MatrixID',
    pattern=NA,
    param_delim='_',
    ...){
    # input_dir=file.path(TAD_DIR, 'method_cooltools'); suffix='-TAD.tsv'; filename.column.name='MatrixID'; pattern=NA; param_delim='_';
    # !!NOTICE!!
    # This will break if any parameter_dir name has a param_delim character in the name or value, not as the delimiter
    # This shouldnt  break if the filename has a single param_delim character in it 
    # i.e. min_resolution-0.45 is fine when param_delim='-'
    #      min_resolution_0.45 is will break with any value for param_delim
    # input_dir=HICREP_DIR; suffix='-hicrep.txt'; filename.column.name='file.pair'; param_delim='_'
    suffix_pattern <- glue('*{suffix}$')
    # List all results files that exist
    input_dir %>% 
    list.files(
        pattern=suffix_pattern,
        recursive=TRUE,
        full.names=FALSE
    ) %>% 
    tibble(fileinfo=.) %>%
    mutate(filepath=file.path(input_dir, fileinfo)) %>% 
    {
        if (is.na(pattern)) {
            .
        } else {
            filter(., grepl(pattern, filepath))
        }
    } %>% 
    # Extract param info from directory names into structured columns
    separate_longer_delim(
        fileinfo,
        delim='/'
    ) %>%
    mutate(
        fileinfo=
            ifelse(
                grepl(suffix_pattern, fileinfo),
                paste(
                    filename.column.name,
                    fileinfo,
                    sep=param_delim
                ) %>% 
                str_remove(suffix),
                fileinfo
            )
    ) %>% 
    separate_wider_delim(
        fileinfo,
        delim='_',
        too_many="merge",  # in case filenames have delim inside
        names=
            c(
                'Parameter',
                'Value'
            )
    ) %>%
    pivot_wider(
        names_from=Parameter,
        values_from=Value
    ) %>% 
    # Fix column types based on content
    readr::type_convert()
}

glue_sample_ID <- function(
    Edit,
    Celltype,
    Genotype,
    CloneID,
    TechRepID,
    prefix='',
    ...) {
    glue('{prefix}{Edit}.{prefix}{Celltype}.{prefix}{Genotype}.{prefix}{CloneID}.{prefix}{TechRepID}')
}

get_info_from_SampleIDs <- function(
    df,
    sample_ID_col='SampleID',
    col_prefix='',
    keep_id=TRUE,
    nest_col=NA,
    SampleID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID'
        ),
    ...){
    df %>%
    mutate(
        !!as.character(glue('{col_prefix}isMerged')) :=
            ifelse(
                grepl('Merged', !!sym(sample_ID_col)),
                'Merged',
                'Individual'
            ) %>% 
            factor()
    ) %>% 
    separate_wider_delim(
        all_of(sample_ID_col),
        delim=fixed('.'),
        names=SampleID.fields %>% paste0(col_prefix, .),
        cols_remove=!keep_id
    ) %>% 
    {
        if (!is.na(nest_col)) {
            if (col_prefix  != '') {
                nest(., !!nest_col := starts_with(col_prefix))
            } else {
                nest(
                    ., 
                    !!nest_col :=
                        all_of(
                            c(
                                SampleID.fields,
                                'isMerged'
                            )
                        )
                )
            }
        } else {
            .
        }
    }
}

get_info_from_MatrixIDs <- function(
    df,
    matrix_ID_col='MatrixID',
    sample_ID_col='SampleID',
    col_prefix='',
    keep_id=TRUE,
    nest_col=NA,
    delim=fixed('.'),
    separate_names=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID',
            NA,
            'ReadFilter',
            NA
        ),
    ...){
    # matrix_ID_col='MatrixID.P1'; sample_ID_col='SampleID.P1'; col_prefix='SampleInfo.P1.'; keep_id=FALSE; nest_col=NA
    # Set up some column names/variables
    col_names <- c(separate_names[!is.na(separate_names)], 'isMerged')
    # extract + format metadata
    df %>%
    mutate(isMerged=grepl('Merged', !!sym(matrix_ID_col))) %>% 
    separate_wider_delim(
        all_of(matrix_ID_col),
        delim=delim,
        names=separate_names,
        too_many='merge',
        cols_remove=!keep_id
    ) %>% 
    # create SampleID
    {
        if (!is.null(sample_ID_col)) {
            mutate( ., !!sample_ID_col := pmap(.l=., .f=glue_sample_ID)) %>%
            unnest(!!sym(sample_ID_col))
        } else {
            .
        }
    } %>% 
    # rename columns with a common prefix
    {
        if (col_prefix != '') {
            rename_with(
                ., 
                .fn=~ paste0(col_prefix, .x), 
                .cols=all_of(col_names)
            )
        } else {
            .
        }
    } %>% 
    # nest and rename columns as specified
    {
        if (!is.na(nest_col) && col_prefix != '') {
            nest(
                .,
                !!nest_col := all_of(paste0(col_prefix, col_names))
            )
        } else if (!is.na(nest_col) && col_prefix == '') {
            nest(
                .,
                !!nest_col := all_of(col_names)
            )
        } else {
            .
        }
    }
}

handle_CLI_args <- function(
    args=c('threads', 'force', 'resolutions'),
    has.positional=FALSE){
    # parse listed arguments
    parsed.args <- 
        OptionParser() %>%
        {
            if ('threads' %in% args){
                add_option(
                    .,
                    c('-t', '--threads'),
                    type='integer',
                    default=length(availableWorkers()),
                    dest='threads'
                )
            } else {
                .
            }
        } %>% 
        {
            if ('force' %in% args){
                add_option(
                    .,
                    c('-f', '--force'),
                    action='store_true',
                    default=FALSE,
                    dest='force.redo'
                )
            } else {
                .
            }
        } %>% 
        {
            if ('resolutions' %in% args){
                add_option(
                    .,
                    c('-r', '--resolutions'),
                    type='character',
                    default=paste(c(10, 25, 50, 100) * 1e3, collapse=','),
                    dest='resolutions'
                )
            } else {
                .
            }
        } %>% 
        parse_args(positional_arguments=TRUE)
    # parse list of resolutions if passed
    if ('resolutions' %in% args){
        parsed.args$options$resolutions <- 
            parsed.args$options$resolutions %>%
            str_split(',') %>%
            lapply(as.integer) %>%
            unlist()
    }
    # return positional args if supplied
    if (has.positional){
        parsed.args
    } else {
        parsed.args$options
    }
}

###################################################
# Format stuff
###################################################
scale_numbers <- function(
    numbers,
    accuracy=2,
    force_numeric=FALSE){
    if ((is.character(numbers) | is.factor(numbers)) | force_numeric) {
        numbers %>%
        as.character() %>% 
        tibble(resolution.str=.) %>%
        mutate(
            suffix=str_extract(resolution.str, '[KMG]b'),
            magnitude=
                case_when(
                    suffix == 'Kb' ~ 1e3,
                    suffix == 'Mb' ~ 1e6,
                    suffix == 'Gb' ~ 1e9,
                    TRUE ~ 1
                ),
            resolution=
                resolution.str %>% 
                str_remove('[KMG]b') %>%
                as.integer() %>%
                multiply_by(magnitude)
        ) %>% 
        pull(resolution) #%>% format(scientific=FALSE) %>% as.numeric()
    } else if (is.numeric(numbers) & !force_numeric) {
        numbers %>%
        tibble(resolution=.) %>%
        mutate(
            magnitude=resolution %>% log10() %>% floor() %>% {. %/% 3} %>% {. * 3},
            suffix=
                case_when(
                    magnitude == 3 ~ 'Kb',
                    magnitude == 6 ~ 'Mb',
                    magnitude == 9 ~ 'Gb',
                    TRUE ~ ''
                ),
            resolution.digits=signif(resolution / 10**magnitude, digits=accuracy),
            resolution.str=glue('{resolution.digits}{suffix}')
        ) %>% 
        mutate(resolution.str=fct_reorder(resolution.str, resolution)) %>% 
        pull(resolution.str)
    } else {
        numbers
    }
}

rename_chrs <- function(
    chrs, 
    to_numeric=FALSE){
    if (is.character(chrs) | to_numeric){
        case_when(
            chrs == 'chrX'     ~ 23,
            chrs == 'chrY'     ~ 24,
            grepl(chrs, 'chr') ~ str_remove(chrs, 'chr'),
            TRUE               ~ NA
        ) %>%
        as.integer()
    } else {
        case_when(
            chrs == 23             ~ 'X',
            chrs == 24             ~ 'Y',
            chrs  >  0 & chrs < 23 ~ as.character(chrs),
            TRUE                   ~ NA
        ) %>% 
        paste0('chr', .) %>% 
        factor(levels=CHROMOSOMES)
    }
}

standardize_data_cols <- function(
    results.df,
    skip.isGenome=FALSE,
    skip.isMerged=FALSE,
    skip.resolution=FALSE,
    skip.window.size=FALSE,
    skip.chr=FALSE){
    results.df %>% 
    {
        if ('isGenome' %in% colnames(.) & !skip.isGenome) {
            mutate(., isGenome=factor(isGenome, levels=c('Per.Chr', 'Genome.Wide')))
        } else {
            .
        }
    } %>% 
    {
        if ('isMerged' %in% colnames(.) & !skip.isMerged) {
            if (is.logical(results.df$isMerged)) {
                mutate(
                    .,
                    isMerged=
                        ifelse(isMerged, 'Merged', 'Individual') %>%
                        factor(levels=c('Merged', 'Individual'))
                )
            } else {
                mutate(., isMerged=factor(isMerged, levels=c('Merged', 'Individual')))
            }
        } else {
            .
        }
    } %>% 
    {
        if ('resolution' %in% colnames(.) & !skip.resolution) {
            mutate(
                .,
                resolution=
                    resolution %>%
                    scale_numbers(force_numeric=TRUE) %>%
                    scale_numbers(),
            )
        } else {
            .
        }
    } %>% 
    {
        if ('window.size' %in% colnames(.) & !skip.window.size) {
            mutate(
                .,
                window.size=
                    window.size %>%
                    scale_numbers(force_numeric=TRUE) %>%
                    scale_numbers(),
            )
        } else {
            .
        }
    } %>% 
    {
        if ('chr' %in% colnames(.) & !skip.chr) {
            if ('Genome.Wide' %in% .$chr) {
                mutate(., chr=factor(chr, levels=c(CHROMOSOMES, 'Genome.Wide'))) %>% 
                mutate(
                    isGenome=
                        ifelse(chr == 'Genome.Wide', 'Genome.Wide', 'Per.Chr') %>% 
                        factor(levels=c('Genome.Wide', 'Per.Chr'))
                )
            } else {
                mutate(., chr=factor(chr, levels=CHROMOSOMES))
            }
        } else {
            .
        }
    }
}

###################################################
# Load Specific Data
###################################################
load_sample_metadata <- function(filter=TRUE){
    SAMPLE_METADATA_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    {
        if(filter) {
            filter(., Included)
        } else {
            .
        }
    }
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=c('Chr', 'chr.total.bp')
    )
}

get_min_resolution_per_matrix <- function(
    df=NULL,
    as_int=TRUE,
    filter_res=TRUE){
    # get minimum viable resolution for each matrix based on Rao et at. 2014 definition
    MIN_SAMPLE_RESOLUTION_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    select(SampleID, resolution) %>% 
    mutate(resolution=scale_numbers(resolution, force_numeric=as_int)) %>% 
    add_column(is.smallest.resolution=TRUE) %>% 
    {
        if (is.null(df)) {
            .
        } else {
            left_join(
                .,
                df,
                by=join_by(SampleID)
            )
        }
    } %>% 
    { 
        if (filter_res) {
            filter(., is.smallest.resolution) %>%
            select(-c(is.smallest.resolution))
        } else {
            .
        }
    }
}

###################################################
# Load Filetypes
###################################################
load_genome_coverage <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
    rename_with(~ str_remove(.x, '_(raw|weight)')) %>% 
    rename(
        'chr'=chrom,
        'coverage.cis'=cov_cis,
        'coverage.total'=cov_tot,
    ) %>% 
    pivot_longer(
        starts_with('coverage'),
        names_to='metric',
        names_prefix='coverage.',
        values_to='coverage'
    ) %>%  
    standardize_data_cols()
}

load_mcool_file <- function(
    filepath,
    resolution=100000,
    normalization="NONE",
    range1="",
    range2="",
    cis=TRUE,
    type='df',
    ...){
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        normalization=normalization,
        type=type,
        join=TRUE,
        query_type='UCSC'
    ) %>% 
    # format column names
    {
        if (type == 'df') {
            if (cis) {
                as_tibble(.) %>%
                filter(chrom1 == chrom2) %>% 
                select(
                    -c(
                        chrom2,
                        end1,
                        end2
                    )
                ) %>% 
                rename(
                    'chr'=chrom1,
                    'range1'=start1,
                    'range2'=start2,
                    'IF'=count
                )
            } else {
                as_tibble(.) %>%
                rename(
                    'chr1'=chrom1,
                    'chr2'=chrom2,
                    'range1'=start1,
                    'range2'=start2,
                    'IF'=count
                )
            }
        } else {
            as.matrix(.)
        }
    }
}

list_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    resolutions=NULL,
    normalizations=NULL,
    ...){
    # List all cooler files 
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    # parse sample metadata
    mutate(MatrixID=basename(filepath)) %>% 
    get_info_from_MatrixIDs(
        matrix_ID_col='MatrixID',
        sample_ID_col='SampleID',
        col_prefix='',
        keep_id=FALSE,
        nest_col=NA,
    ) %>% 
    # List all paramter combinations
    {
        if (!is.null(normalizations)) {
            cross_join(., tibble(normalization=normalizations))
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(resolutions)) {
            cross_join(., tibble(resolution=resolutions))
        } else {
            .
        }
    }
}

load_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    resolutions=NULL,
    normalizations=NULL,
    regions.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    return_metadata_only=FALSE,
    keep_metadata_columns=FALSE,
    ...){
    # Define all genomic regions to load contacts 
    regions.df <- 
        if (is.null(regions.df)) {
            # Get the whole genome
            if ((is.null(range1s)) & (is.null(range2s))) {
                tibble()
                # tibble(
                #     range1=CHROMOSOMES,
                #     range2=CHROMOSOMES
                # )
            # get all contacts within all regions in range1s
            } else if (is.null(range2s)) {
                tibble(
                    range1=range1s,
                    range2=range1s
                )
            # get all contacts between all pairs of regions only (not intra-region contacts)
            } else {
                expand_grid(
                    range1=range1s,
                    range2=range2s
                )
            }
        # Just load the specified regions (1 region per row: chr, start, end)
        } else {
            regions.df
        }
    # List all regions for all samples
    list_mcool_files(
        pattern=pattern,
        resolutions=resolutions,
        normalizations=normalizations,
    ) %>% 
    join_all_rows(regions.df) %>% 
    # Load contacts if specified or just return sample metadata + filepaths + regions
    {
        if (return_metadata_only) {
            .
        } else {
            mutate(
                .,
                contacts=
                    purrr::pmap(
                        .l=.,
                        .f=load_mcool_file,
                        .progress=progress
                    )
            ) %>% 
            select(-c(filepath, range1, range2)) %>% 
            unnest(contacts)
        }
    }
}

###################################################
# Group samples by condition
###################################################
set_up_sample_groups <- function(
    sample.groups,
    resolutions=c(),
    use_merged=FALSE){
    list_mcool_files() %>%
    get_min_resolution_per_matrix() %>% 
    distinct() %>% 
    {
        if (use_merged) {
            filter(., isMerged)
        } else {
            filter(., !isMerged)
        }
    } %>% 
    # Now group samples by condition, 
    nest(samples.df=-c(isMerged)) %>% 
    cross_join(
        sample.groups %>% 
        mutate(
            across(
                starts_with('Sample.Group'),
                ~ str_replace_all(.x, 'All', '.*'),
                .names='{.col}.Pattern'
            )
        )
    ) %>%
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.Pattern) ~ Sample.Group,
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
    ungroup() %>%
    mutate(resolution=list(unique(c(resolutions, resolution.min, resolution.max)))) %>%
    unnest(resolution) %>% 
    mutate(
        resolution.type=
            case_when(
                resolution == resolution.max ~ 'max',
                resolution == resolution.min ~ 'min',
                TRUE                         ~ NA
            )
    ) %>% 
    select(
        -c(
            resolution.min,
            resolution.max,
            ends_with('.Pattern')
        )
    )
    
}

set_up_sample_comparisons <- function(
    comparison.groups,
    resolutions=c(),
    merging='individual'){
    # get info + filepaths for all contact matrices
    list_mcool_files() %>%
    get_min_resolution_per_matrix() %>% 
    distinct() %>% 
    {
        if (merging == 'individual' ) {
            filter(., !isMerged)
        } else if (merging == 'merged' ) {
            filter(., isMerged)
        } else {
            .
        }
    } %>% 
    # Now group samples by condition, 
    nest(samples.df=-c(isMerged)) %>% 
    cross_join(
        comparison.groups %>% 
        mutate(
            across(
                starts_with('Sample.Group.'),
                ~ str_replace_all(.x, 'All', '.*'),
                .names='{.col}.Pattern'
            )
        )
    ) %>%
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.Left.Pattern) ~ Sample.Group.Left,
                        str_detect(SampleID, Sample.Group.Right.Pattern) ~ Sample.Group.Right,
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
    # List every comparison + resolution that is either a min or max for 1 comparison
    ungroup() %>%
    mutate(resolution=list(unique(c(resolutions, resolution.min, resolution.max)))) %>%
    unnest(resolution) %>% 
    mutate(
        resolution.type=
            case_when(
                resolution == resolution.max ~ 'max',
                resolution == resolution.min ~ 'min',
                TRUE                         ~ NA
            )
    ) %>% 
    select(-c(resolution.min, resolution.max)) %>%
    select(-c(ends_with('.Pattern')))
}

sample_group_priority_fnc_16p <- function(Sample.Group){
    # FC is determined by the edger::exactTest() function called in multiHiCCompare
    # https://github.com/dozmorovlab/multiHiCcompare/blob/dcfe4aaa8eaef45e203f3d7f806232bb613d2c9b/R/glm.R#L69
    # According to the docs for exactTest()
    # """Note that the first group listed in the pair is the baseline for the comparison—so if the pair is c("A","B") then the comparison is B - A, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change)."""
    # So with the factor level that comes FIRST is used as the baseline i.e. DENOMINATOR
    # https://www.quantargo.com/help/r/latest/packages/edgeR/NEWS/exactTest
    # so for 16p.NSC.DEL vs 16p.NSC.WT we want to force 16p.NSC.DEL be numerator 
    # therefore we make it the SECOND factor level i.e. have  a larger priority number 
    # i.e. this works when the priority of 16p.NSC.DEL > 16p.NSC.WT
    # so if we create a factor based on this priority the factor hass the following levels
    # Levels: 16p.NSC.WT 16p.NSC.DEL
    # which results in the call exactTest(groups=c(16p.NSC.WT, 16p.NSC.DEL))
    # which is DEL as the FC numerator
    case_when(
        Sample.Group == '16p.NSC.DUP' ~ 1,  # always numerator in FCs
        Sample.Group == '16p.iN.DUP'  ~ 2,
        Sample.Group == '16p.NSC.DEL' ~ 3,
        Sample.Group == '16p.iN.DEL'  ~ 4,
        Sample.Group == '16p.NSC.WT'  ~ 5,
        Sample.Group == '16p.iN.WT'   ~ 6,
        TRUE                          ~ Inf
    )
}

sample_group_priority_fnc_Cohesin <- function(Sample.Group){
    # FC is determined by the edger::exactTest() function called in multiHiCCompare
    # https://github.com/dozmorovlab/multiHiCcompare/blob/dcfe4aaa8eaef45e203f3d7f806232bb613d2c9b/R/glm.R#L69
    # According to the docs for exactTest()
    # """Note that the first group listed in the pair is the baseline for the comparison—so if the pair is c("A","B") then the comparison is B - A, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change)."""
    # So with the factor level that comes FIRST is used as the baseline i.e. DENOMINATOR
    # https://www.quantargo.com/help/r/latest/packages/edgeR/NEWS/exactTest
    # so for NIPBL.DEL vs NIPBLWT we want to force NIPBL.DEL to be numerator therefore we make it 
    # the SECOND factor level i.e. have  a larger priority number 
    # i.e. this works when the priority of NIPBL.DEL > NIPBL.WT 

    case_when(
        grepl( 'CTCF.iN.DEL', Sample.Group) ~ 1, # always numerator since always last factor level
        grepl('RAD21.iN.DEL', Sample.Group) ~ 2,
        grepl( 'WAPL.iN.DEL', Sample.Group) ~ 3,
        grepl('NIPBL.iN.DEL', Sample.Group) ~ 4,
        grepl(  'All.iN.DEL', Sample.Group) ~ 5, 
        grepl( 'CTCF.iN.WT',  Sample.Group) ~ 11,
        grepl('RAD21.iN.WT',  Sample.Group) ~ 12,
        grepl( 'WAPL.iN.WT',  Sample.Group) ~ 13,
        grepl('NIPBL.iN.WT',  Sample.Group) ~ 14,
        grepl(  'All.iN.WT',  Sample.Group) ~ 15,
        TRUE                                ~ Inf
    )
}

set_foldchange_direction_as_factor <- function(
    results.df,
    sample_group_priority_fnc,
    group1_colname='Sample.Group.Left', 
    group2_colname='Sample.Group.Right',
    ...){
    # Use this to make sure that when testing is done by edger::exactTest(), which relies only 
    # on factor levels to determine fold-change direction, that the explicitly stated numerator
    group1_priority <- glue('{group1_colname}.Priority') 
    group2_priority <- glue('{group2_colname}.Priority') 
    # will have the approriate factor level detected by edger
    results.df %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    # use priority function to determine which sample group should represent the numerator 
    # in the fold change 
    mutate(
        across(
            c(group1_colname, group2_colname),
            ~ sample_group_priority_fnc(.x),
            .names='{.col}.Priority'
        )
    ) %>%
    # Create explicit and consistent Numerator column, i.e. if there is a log(FC) > 0 then 
    # the Numerator is enriched for the signal being compared
    # The Denominator is set to be the opposite group
    mutate(
        Sample.Group.Numerator=
            case_when(
                .data[[group1_priority]] < .data[[group2_priority]] ~ .data[[group1_colname]],
                .data[[group1_priority]] > .data[[group2_priority]] ~ .data[[group2_colname]],
                TRUE ~ NA
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == .data[[group1_colname]] ~ .data[[group2_colname]],
                Sample.Group.Numerator == .data[[group2_colname]] ~ .data[[group1_colname]],
                TRUE ~ NA
            )
    ) %>% 
    # clean up unneded columns since we now have explicity numerator/denominator labels
    select(
        -c(
            ends_with('.Priority'),
            group1_colname,
            group2_colname
        )
    )
}

###################################################
# Handle pairs of samples
###################################################
join_all_rows <- function(
    df1,
    df2=NULL,
    cols_to_match=c(),
    ...){
    if (is.null(df2)){ 
        df2 <- df1
    }
    if (length(cols_to_match) == 0) {
        cross_join(
            df1,
            df2,
            ...
        ) 
    } else {
        full_join(
            df1,
            df2,
            relationship='many-to-many',
            by=cols_to_match,
            ...
        )
    }
}

get_all_row_combinations <- function(
    df1,
    df2=NULL,
    col_to_pair=NA,
    cols_to_match=c(),
    suffixes=c('.P1', '.P2'),
    keep_self=TRUE,
    ...){
    # df2=NULL; col_to_pair='SampleID'; cols_to_match=c('resolution', 'chr', 'TAD.method', 'TAD.params'); suffixes=c('.P1', '.P2'); keep_self=FALSE;
    # Get all combinations of rows with matching attributes (cols_to_pair)
    colA <- glue('{col_to_pair}.A')
    colB <- glue('{col_to_pair}.B')
    join_all_rows(
        df1,
        df2,
        cols_to_match=cols_to_match,
        suffix=c('.A', '.B')
    ) %>% 
    # keep/remove pairs of a row matched to itself
    {
        if (keep_self) {
            .
        } else {
            filter(., !!sym(colA) != !!sym(colB))
        }
    } %>% 
    # Keep only distinct pairs of rows (order doesnt matter)
    rowwise() %>% 
    mutate(
        index=
            case_when(
                !!sym(colA) >= !!sym(colB) ~ paste(!!sym(colA), !!sym(colB), collapse='~'),
                !!sym(colA) <  !!sym(colB) ~ paste(!!sym(colB), !!sym(colA), collapse='~')
            ) %>% 
            paste(collapse='~')
    ) %>% 
    ungroup() %>% 
    group_by(across(c(cols_to_match))) %>% 
    distinct(
        index,
        .keep_all=TRUE
    ) %>%
    ungroup() %>% 
    select(-c(starts_with('index'))) %>%
    # rename with suffixes
    rename_with(~ str_replace(.x, '\\.A$', suffixes[[1]])) %>% 
    rename_with(~ str_replace(.x, '\\.B$', suffixes[[2]]))
}

merge_sample_info <- function(
    SampleInfo.P1,
    SampleInfo.P2,
    suffix.P1='.P1',
    suffix.P2='.P2',
    ...){
    # SampleInfo.P1=df$SampleInfo.P1[[1]]; prefix.P1='SampleInfo.P1.'; SampleInfo.P2=df$SampleInfo.P2[[1]]; prefix.P2='SampleInfo.P2.'
    # Combine pair info
    if (is.null(SampleInfo.P1)) {
        return(NULL)
    } else if (is.null(SampleInfo.P2)) {
        return(NULL)
    }
    bind_rows(
        SampleInfo.P1 %>% rename_with(~ str_remove(.x, suffix.P1)),
        SampleInfo.P2 %>% rename_with(~ str_remove(.x, suffix.P2)),
    ) %>% 
    add_column(PairIndex=c('P1', 'P2')) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(
        -PairIndex,
        names_to='SampleAttribute',
        values_to='SampleValue'
    ) %>% 
    pivot_wider(
        names_from=PairIndex,
        values_from=SampleValue
    ) %>% 
    rowwise() %>% 
    # Now store both pair values for a single metadata field in 1 column
    # and keep order consistent for grouping downstream
    # NOTE this means that the metadata columns are not always P1 vs P2
    # sometimes it will be P2 vs P1 based on alphabetical order of the values
    # but thats fine since you can just looks at the SampleIDs themselves if 
    # you care about that info
    mutate(
        PairValue=
            case_when(
                P1 == P2 ~ glue('{P1} vs {P2}'),
                P1 != P2 ~ 
                    c(P1, P2) %>% 
                    sort() %>%  
                    paste(collapse=" vs ")
            )
    ) %>%
    select(SampleAttribute, PairValue) %>%
    pivot_wider(
        names_from=SampleAttribute,
        values_from=PairValue
    )
}

