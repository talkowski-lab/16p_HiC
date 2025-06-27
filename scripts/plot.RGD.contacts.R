##################
## Dependencies
##################
library(here)
here::i_am('scripts/plot.RGD.contacts.R')
BASE_DIR <- here()
# BASE_DIR <- here('../remote.wapl')
SCRIPT_DIR <- file.path(BASE_DIR, 'scripts')
source(file.path(SCRIPT_DIR, 'locations.R'))
# SCRIPT_DIR <- "~/TalkowskiLab/Projects/HiC/remote.16p/scripts"
source(file.path(SCRIPT_DIR, 'constants.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
library(tidyverse)
library(magrittr)
library(HiContacts)
library(furrr)
# PLOT_DIR <- file.path(RESULTS_DIR, 'plots/RGDs')
PLOT_DIR <- file.path(RESULTS_DIR, 'plots/lab.meeting')
##################
# Load all fileinfo + regions, but not contact data itself
##################
# RGD regions
# Genomic deletion regions
annotated.contacts.df <- 
    LAB_MEETING_REGIONS %>% 
    fetch_regions(
        normalization=
            c(
                "NONE",
                "weight"
            ),
        resolutions=
            c(
                100e3,
                50e3,
                25e3,
                10e3,
                5e3
            )
    ) %>% 
    format_plot_params() #%>% 
    # get_min_resolution_per_matrix(
    #     int_res=TRUE,
    #     filter_res=TRUE
    # )
# plot shape
plot.params <- 
    tribble(
~plot_type,         ~make_diagonal,  ~make_symmetric,  ~add_NAs, ~ylinewidth, ~height,
'square.heatmap',            FALSE,             TRUE,      TRUE,         0.7,       8,
# 'diagonal.heatmap',           TRUE,             TRUE,     FALSE,           0,       8,
# 'triangle.heatmap',           TRUE,            FALSE,     FALSE,           0,       6
)
# annotated.contacts.df %>% count(Sample.ID, resolution)
# RGD regions
# annotated.contacts.df <- 
#     LAB_MEETING_REGIONS %>% 
#     fetch_regions(
#         normalization=
#             c(
#                 "NONE",
#                 "weight"
#             ),
#         resolutions=
#             c(
#                 50e3,
#                 5e3
#             )
#     ) %>% 
#     format_plot_params()
##################
# Plot contact heatmaps
##################
annotated.contacts.df %>% 
    # filter(region == '17q21.31') %>% 
    join_all_rows(plot.params) %>% 
    filter(region == '16p11.2 CNV') %>%  
    filter(grepl('NSC', Sample.ID)) %>% 
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
                plot_type,
                glue('{Sample.ID}-{plot_type}.pdf')
            )
    ) %>% 
    # 1 plot (file) per param combo
    pmap(
        .l=.,
        .f=heatmap_wrapper,
        x_var='range1',
        y_var='range2',
        fill_var='IF',
        transform_fnc=log10,
        cmap=coolerColors(),
        na.color='white',
        width=10,
        .progress=TRUE
    )
##################
# Plot contact heatmaps and log2(A/B) heatmaps for regions of interest
##################
annotated.contacts.df %>%
    mutate(isMerged=grepl('Merged', Sample.ID)) %>% 
    make_sample_pairs() %>% 
    mutate(
        make_symmetric=TRUE,
        add_explicity_empty_bins=TRUE,
        x_var='range1',
        y_var='range2',
        fill_var='IF',
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
                fill_var,
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
            fill_var,
            range1,
            range2,
            cis
        )
    ) %>% 
    pmap(
        .l=.,
        .f=logfc_heatmap_wrapper,
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
# plot entire chromosomes
    pmap(
        .l=.,
        .f=heatmap_wrapper,
        x_var='range1',
        y_var='range2',
        fill_var='IF',
        transform_fnc=log10,
        cmap=coolerColors(),
        na.color='white',
        width=10,
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
