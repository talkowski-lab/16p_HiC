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
# Load all fileinfo + function params, but not contact data itself
##################
annotated.contacts.df <- 
    # regions of interest to plot
    fetch_RGD_regions(
        normalization=
            c(
                "NONE",
                "weight"
            ),
        resolutions=
            c(
              5e4,
              5e3
            )
    ) %>% 
    # organize all files + params + coords, 1 plot per file
    format_rgd_plot_params()
##################
# Plot square contact heatmap
##################
annotated.contacts.df %>% 
    # filter(normalization == 'NONE', grepl('WAPL.DEL.Merged', Sample.ID)) %>% 
    # Set plot labels/numbers + output_file etc.
    mutate(
        xlab=glue('{region.chr} Position'),
        ylab=glue('{region.chr} Position'),
        fill_lab=
            case_when(
                normalization == 'NONE' ~ 'log10(contacts)',
                normalization == 'weight' ~ 'log10(balanced)',
            ),
        output_file=
            file.path(
                output_dir,
                'square.heatmap', 
                glue('{Sample.ID}-square.heatmap.pdf')
            )
    ) %>% 
    # # 1 plot (file) per param combo
    group_by(output_file) %>% 
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
         na.value="grey50",
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
    filter(normalization == 'NONE', grepl('WAPL.DEL.Merged', Sample.ID)) %>% 
    # Set plot labels/numbers + output_file etc.
    mutate(
        xlab=glue('{region.chr} Position'),
        ylab=glue('{region.chr} Position'),
        fill_lab=
            case_when(
                normalization == 'NONE' ~ 'log10(contacts)',
                normalization == 'weight' ~ 'log10(balanced)',
            ),
        output_file=
            file.path(
                output_dir,
                'triangle.heatmap', 
                glue('{Sample.ID}-triangle.heatmap.pdf')
            )
    ) %>% 
    # # 1 plot (file) per param combo
    group_by(output_file) %>% 
    pmap(
         .l=.,
         .f=heatmap_wrapper,
         x_var='range1',
         y_var='range2',
         fill_var='IF',
         transform_fnc=log10,
         cmap=coolerColors(),
         xlinewidth=0.7,
         ylinewidth=0.0,
         width=10,
         height=8,
         na.value="grey50",
         .progress=TRUE
    )
##################
# TESTING STUFF
##################
# annotated.contacts.df  %>% distinct(window.size, resolution, range1, region)
# annotated.contacts.df$output_dir[1:5]
# annotated.contacts.df %>% colnames()
# annotated.contacts.df %>% head(2) %>% t()
# annotated.contacts.df %>% count(Sample.ID, region, normalization)
# sample.pair.df %>% colnames()
# sample.pair.df$contacts[[1]] %>% colnames()
