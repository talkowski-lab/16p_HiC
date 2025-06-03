##################
## Dependencies
##################
library(here)
here::i_am('scripts/plot.RGD.contacts.R')
BASE_DIR <- here()
SCRIPT_DIR <- file.path(BASE_DIR, 'scripts')
# SCRIPT_DIR <- "~/TalkowskiLab/Projects/HiC/remote.16p/scripts"
source(file.path(SCRIPT_DIR, 'constants.R'))
source(file.path(SCRIPT_DIR, 'locations.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
library(tidyverse)
library(magrittr)
library(HiContacts)
library(furrr)
PLOT_DIR <- file.path(RESULTS_DIR, 'plots/RGDs')
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
# Plot contact heatmaps and log2(A/B) heatmaps for regions of interest
##################
# plan(multisession, workers=16)
annotated.contacts.df %>%
    mutate(isMerged=grepl('Merged', Sample.ID)) %>% 
    make_sample_pairs() %>% 
    mutate(
        make_symmetric=TRUE,
        add_explicity_empty_bins=FALSE,
        x_var='range1',
        y_var='range2',
        xlab=glue('{region.chr} Position'),
        ylab=glue('{region.chr} Position'),
        output_file=
            file.path(
                output_dir,
                'logfc.heatmap',
                glue('{Sample.ID.A}_vs_{Sample.ID.B}-logfc.heatmap.pdf')
            )
    ) %>% 
    nest(
        contact.params=
            c(
                make_symmetric,
                add_explicity_empty_bins,
                x_var,
                y_var,
                resolution,
                normalization,
                range1,
                range2,
                cis
            )
    ) %>% 
    mutate(wide_params=contact.params) %>%
    unnest(wide_params) %>% 
    select(
        -c(
            make_symmetric,
            add_explicity_empty_bins,
            x_var,
            y_var,
            range1,
            range2,
            cis
        )
    ) %>% 
    group_by(output_file) %>% 
    pmap(
        .l=.,
        .f=logfc_heatmap_wrapper,
        fill_var='log2fc',
        fill_lab="log2(A / B)",
        sample_priority_fnc=NIPBL_WAPL_sample_priority_fnc,
        transform_fnc=function(x) 1 * x,
        cmap=bgrColors(),
        na.value='grey50',
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
