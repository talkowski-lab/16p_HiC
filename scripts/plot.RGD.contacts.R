##################
## Dirs
##################
library(here)
here::i_am('scripts/plot.RGD.contacts.R')
BASE_DIR=here()
source(file.path(BASE_DIR, 'scripts/locations.R'))
# SCRIPT_DIR <- "/home/sid/TalkowskiLab/Projects/HiC/remote.16p/scripts"
PLOT_DIR <- file.path(RESULTS_DIR, 'plots/RGDs')
##################
## Dependencies
##################
library(tidyverse)
library(magrittr)
library(HiContacts)
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
##################
# regions of interest to plot
##################
region.df <- 
    load_rgd_regions(
        normalization=c("NONE", "weight"),
        window.size=1e6
    )
##################
# Load all contacts within eac hregion
##################
annotated.contacts.df <- 
    check_cached_results( 
        results_file=RGD_ANNOTATED_CONTACTS_FILE,
        force_redo=TRUE,
        return_data=TRUE,
        results_fnc=load_rgd_contacts,
        region.df=region.df,
        resolutions=c(5e4)
    ) %>% 
    mutate(
        normalization=
            case_when(
                normalization == 'weight' ~ 'ICE',
                TRUE ~ normalization
            ),
        region.title=glue('{region} RGD Region {region.UCSC}'),
        across(
            c(resolution, window.size),
            scale_numbers
        ),
        output_dir=
            file.path(
                PLOT_DIR, 
                glue('region_{region}'),
                glue('normalization_{normalization}'),
                glue('resolution_{resolution}'),
                glue('context_{window.size}')
            )
    )
##################
# Plot square contact heatmap
##################
annotated.contacts.df %>% 
    # Make symmetrical 
    bind_rows(
        annotated.contacts.df %>% 
        rename('range2'=range1, 'range1'=range2) %>%
        filter(range1 != range2)
    ) %>% 
    # Make empty bin-pairs explicit NAs
    group_by(across(-c(range1, range2, IF))) %>% 
    complete(range1, range2, fill=list(IF=NA)) %>% 
    ungroup() %>%
    nest(contacts=c(region.start, region.end, range1, range2, IF)) %>% 
    mutate(
        xlab=glue('{chr} Position'),
        ylab=glue('{chr} Position'),
        fill_lab=
            case_when(
                normalization == 'NONE' ~ 'log10(contacts)',
                normalization == 'ICE' ~ 'log10(balanced)',
            ),
        output_file=
            file.path(
                output_dir,
                'square.heatmap', 
                glue('{Sample.ID}-square.heatmap.pdf')
            )
    ) %>% 
    group_by(Sample.ID, region.title) %>% 
    pmap(
         .l=.,
         .f=heatmap_wrapper,
         x_var='range1',
         y_var='range2',
         fill_var='IF',
         transform_fnc=log10,
         cmap=coolerColors(),
         axis_label_accuracy=0.01,
         x_text_angle=0,
         xlinewidth=0.7,
         ylinewidth=0.7,
         width=10,
         height=8,
         .progress=TRUE
    )
##################
# Get all pairs of samples for all possible plots
##################
sample.pair.df <- 
    annotated.contacts.df %>%
    bind_rows(
        annotated.contacts.df %>% 
        rename('range2'=range1, 'range1'=range2) %>%
        filter(range1 != range2)
    ) %>% 
    mutate(
        isMerged=grepl('Merged', Sample.ID),
        region.chr=chr
    ) %>% 
    # region specific info
    nest(
        region.info=
            c(
                output_dir,
                region.title,
                region.chr,
                region.UCSC,
                region.dist
            )
    ) %>% 
    # Sample+region specific info
    nest(
        contacts=
            c(
                region.start,
                region.end,
                chr,
                range1,
                range2,
                IF
            )
    ) %>% 
    # Get all pairs of samples with comparable contact matrices
    full_join(
        .,
        {.},
        suffix=c('.A', '.B'),
        copy=TRUE,
        by= # things that must be consistent between samples being compared
            join_by(
                isMerged,
                ReadFilter,
                region,
                normalization,
                resolution,
                window.size
            )
    ) %>%
    # Dedup pairs (symetrical)
    rowwise() %>% 
    mutate(tmp=paste0(sort(c(Sample.ID.A, Sample.ID.B)), collapse="")) %>% 
    ungroup() %>% 
    distinct(
        tmp, 
        region, 
        ReadFilter,
        normalization,
        resolution,
        window.size,
        .keep_all=TRUE
    ) %>%
    filter(Sample.ID.A != Sample.ID.B) %>% 
    # Also make sure Sample.ID.A is always the Deletion or WAPL sample (FC numerator)
    mutate(
        Sample.ID.numerator=
            case_when(
                grepl('WAPL.DEL', Sample.ID.A)  ~ Sample.ID.A,
                grepl('WAPL.DEL', Sample.ID.B)  ~ Sample.ID.B,
                grepl('NIPBL.DEL', Sample.ID.A) ~ Sample.ID.A,
                grepl('NIPBL.DEL', Sample.ID.B) ~ Sample.ID.B,
                grepl('WAPL.WT', Sample.ID.A)   ~ Sample.ID.A,
                grepl('WAPL.WT', Sample.ID.B)   ~ Sample.ID.B,
                TRUE ~ Sample.ID.A
            ),
        Sample.ID.denominator=
            case_when(
                Sample.ID.numerator == Sample.ID.A ~ Sample.ID.B,
                Sample.ID.numerator == Sample.ID.B ~ Sample.ID.A
            ),
        contacts.numerator=
            case_when(
                Sample.ID.numerator == Sample.ID.A ~ contacts.A,
                Sample.ID.numerator == Sample.ID.B ~ contacts.B
            ),
        contacts.denominator=
            case_when(
                Sample.ID.numerator == Sample.ID.A ~ contacts.B,
                Sample.ID.numerator == Sample.ID.B ~ contacts.A
            )
    ) %>% 
    select(
        -c(
            Sample.ID.A,
            Sample.ID.B,
            contacts.A,
            contacts.B
        )
    ) %>% 
    rename(
        'Sample.ID.A'=Sample.ID.numerator,
        'Sample.ID.B'=Sample.ID.denominator,
        'contacts.A'=contacts.numerator,
        'contacts.B'=contacts.denominator
    ) %>% 
    # make plot df for each pair
    rowwise() %>% 
    mutate(
        contacts=
            inner_join(
                contacts.A,
                contacts.B,
                suffix=c('.A', '.B'),
                by=
                    join_by(
                        region.start,
                        region.end,
                        chr,
                        range1,
                        range2
                    )
            ) %>%
            mutate(log2fc=log2(IF.A / IF.B)) %>% 
            # group_by(across(-c(range1, range2, IF.A, IF.B, log2fc))) %>% 
            # complete(range1, range2, fill=list(IF=NA)) %>% 
            # ungroup() %>%
            list()
    ) %>% 
    ungroup() %>% 
    # clean up
    unnest(region.info.A) %>% 
    select(
        -c(
            region.info.B,
            tmp,
            contacts.A,
            contacts.B
        )
    )
##################
# Plot contact heatmaps and log2(A/B) heatmaps for regions of interest
##################
sample.pair.df %>%
    # filter(grepl(Sample.ID, 'WAPL.P2C4')) %>% 
    rename('chr'=region.chr) %>% 
    mutate(
        Sample.ID=glue('{Sample.ID.A} vs {Sample.ID.B}'),
        xlab=glue('{chr} Position'),
        ylab=glue('{chr} Position'),
        output_file=
            file.path(
                output_dir,
                'logfc.heatmap',
                glue('{Sample.ID.A}_vs_{Sample.ID.B}-logfc.heatmap.pdf')
            )
    ) %>% 
    group_by(Sample.ID, region.title) %>% 
    pmap(
         .l=.,
         .f=heatmap_wrapper,
         x_var='range1',
         y_var='range2',
         fill_var='log2fc',
         fill_lab="log2(A / B)",
         transform_fnc=function(x) 1 * x,
         cmap=bgrColors(),
         axis_label_accuracy=0.01,
         x_text_angle=25,
         xlinewidth=0.5,
         ylinewidth=0.5,
         width=10,
         height=8,
         .progress=TRUE
    )
##################
# Transform to position vs distance to plot rectangle matrix
##################
annotated.contacts.df %>%
    filter(grepl('WAPL.DEL.Merged', Sample.ID)) %>% 
    # Make empty bin-pairs explicit NAs
    bind_rows(
        annotated.contacts.df %>% 
        rename('range2'=range1, 'range1'=range2) %>%
        filter(range1 != range2)
    ) %>% 
    # Transoforming for triangle version, stolen from HiContacts library
    # https://github.com/js2264/HiContacts/blob/146485f3e517429f7d4aa408898f6e89d4a04b36/R/plotting.R#L422
    mutate(
        resolution.int=scale_numbers(resolution, toint=TRUE),
        window.size.int=scale_numbers(window.size, toint=TRUE),
        # center1=range1 + (resolution.int / 2),
        # center2=range2 + (resolution.int / 2),
        # distance=abs(center2 - center1),
        distance=abs(range2 - range1),
        # coord=range1
        coord=range1 + (range2 - range1) / 2,
        # coord=center1 + (center2 - center1) / 2,
    ) %>% 
    select(-c(range1, range2)) %>% 
    nest(contacts=c(region.start, region.end, coord, distance, IF)) %>% 
    mutate(
        xlab=glue('{chr} Position'),
        ylab=glue('Distance'),
        fill_lab=
            case_when(
                normalization == 'NONE' ~ 'log10(contacts)',
                normalization == 'ICE' ~ 'log10(balanced)',
            ),
        output_file=
            file.path(
                output_dir,
                'distance.heatmap', 
                glue('{Sample.ID}-distance.heatmap.pdf')
            )
    ) %>% 
    group_by(Sample.ID, region.title) %>% 
    pmap(
         .l=.,
         .f=heatmap_wrapper,
         x_var='coord',
         y_var='distance',
         fill_var='IF',
         transform_fnc=log10,
         cmap=coolerColors(),
         axis_label_accuracy=0.01,
         x_text_angle=25,
         na.value="white",
         xlinewidth=0.7,
         ylinewidth=0,
         cfr=(1/sqrt(2)/sqrt(2)),
         width=10,
         height=7,
         .progress=TRUE
    )
##################
# TESTING STUFF
##################
# annotated.contacts.df 
# annotated.contacts.df$output_dir[1:5]
# annotated.contacts.df %>% colnames()
# annotated.contacts.df %>% head(2) %>% t()
# annotated.contacts.df %>% count(Sample.ID, region, normalization)
# sample.pair.df %>% colnames()
# sample.pair.df$contacts[[1]] %>% colnames()
