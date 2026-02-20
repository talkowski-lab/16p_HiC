---
title: "Replicate Reuslts from Weiner et al. 2022"
author: "Siddharth Reed"
date: "`r Sys.Date()`"
output: 
    html_document:
        keep_md: yes
        toc: true
        toc_float: true
        theme: paper
        code_folding: hide
        df_print: paged
---
<style type="text/css">
.main-container {
  max-width: 1800px !important;
  margin-left: auto;
  margin-right: auto;
}
</style>

# Set Up

```{r knitr}
library(knitr)
knitr::opts_chunk$set(
    dev=c('png', 'pdf', 'tiff'),
    dpi=300
)
```
Set up file locations + load dependenies
```{r dependencies, results=FALSE, message=FALSE, warning=FALSE}
library(here)
here::i_am('notebooks/replicate.robinson.results.Rmd')
BASE_DIR <- here()
SCRIPT_DIR <- here('scripts')
source(file.path(SCRIPT_DIR, 'locations.R'))
source(file.path(SCRIPT_DIR, 'constants.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
source(file.path(SCRIPT_DIR, 'utils.robinson.R'))
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(purrr)
```

# Specify regions of interest

We want to know if a specific pair of genomic regions has more contacts than expected compared to randomly selected regions pairs.
Specifically we want to comapre if contacts between the chr16p CNV region and the chr16 telomere are statistically more frequent than between any other regions on chr16p.
```{r regions_of_interest}
# regions to annotate and compare
regions.of.interest <- 
    GENOMIC_REGIONS %>%
    filter(
        region %in% c(
            "chr16p.deletion",
            "chr16p.telomere",
            "chr16p",
            "chr16q"
        )
    ) %>% 
    mutate(region.dist=scale_numbers(region.dist))
# Print regions
regions.of.interest %>%
    select(c(region, region.chr, region.start, region.end, region.dist, region.UCSC)) %>% 
    rename(
        'Region'=region,
        'Chr'=region.chr,
        'Start Pos'=region.start,
        'End Pos'=region.end,
        'Size'=region.dist,
        'UCSC'=region.UCSC
    ) %>% 
    knitr::kable()
```
We also only compare distance-matched contacts to the CNV-telomere contacts of interest, using the following distance range of those CNV-telomere contacts
```{r distance_range}
# Boundaries of regions to compare
telomere.start <- 
    {regions.of.interest %>% filter(region == 'chr16p.telomere') %>% pull(region.start)}
telomere.end <- 
    {regions.of.interest %>% filter(region == 'chr16p.telomere') %>% pull(region.end)}
deletion.start <- 
    {regions.of.interest %>% filter(region == 'chr16p.deletion') %>% pull(region.start)}
deletion.end <- 
    {regions.of.interest %>% filter(region == 'chr16p.deletion') %>% pull(region.end)}
# Distance ranges to compare
dist.min <- 
    round(
        deletion.start - telomere.end,
        digits=-5 # round to nearest 100Kb bin
    )
dist.max <- 
    round(
        deletion.end - telomere.start,
        digits=-5 # round to nearest 100Kb bin
    )
# print
glue("Telomere region is between {telomere.start}-{telomere.end}")
glue("Deletion region is between {deletion.start}-{deletion.end}")
glue("All control contacts used are between {dist.min}bp and {dist.max}bp apart") %>% message()
```

# Load + Annotate Contact Data

Load contact matrices for each file 
```{r load_contacts, rows.print=5}
# Load all contacts for all samples + annotated contacts for which region they are in
contacts.df <- 
    check_cached_results(
        results_file=ROBINSON_REPLICATION_DATA_FILE,
        # force_redo=TRUE,
        results_fnc=load_mcool_files,
        pattern='.hg38.mapq_30.1000.mcool',
        range1s='chr16',
        # resolutions=c(25, 50, 100) * 1e3,
        resolutions=c(25, 50) * 1e3,
        normalizations="NONE",
        cis=TRUE
    ) %>% 
    mutate(resolution=scale_numbers(resolution)) %>% 
    select(-c(ReadFilter))
contacts.df
contacts.df %>% count(resolution, normalization, SampleID, chr)
```
Annotate region-pair of interest for each bin-pair (contacts, 1 per row)
```{r annotate_contacts, rows.print=5}
annotated.contacts.df <- 
    contacts.df %>% 
    annotate_contact_region_pairs(
        regions.of.interest=regions.of.interest,
        most_specific_only=TRUE,
    )
```
Now format annotated contact data for easier plotting
```{r format_annotated_contacts_for_plotting}
# Keep only distance matched control bin-pairs and specify names+colors
plot.df <- 
    annotated.contacts.df %>% 
    format_annotated_contacts()
plot.df %>% dplyr::count(IF.Threshold, Sample.ID, region.title)
```

# Plot IFs of Deleted CNV vs Controls regions 

Now we prepare the data to replicate the results in [Fig.4C from Weiner et al. 2022](https://www.nature.com/articles/s41588-022-01203-y#Fig4).
We use the same procedure to define "Distance-matched control" regions to compare to our region of interest.
Specifically, controls are all pairs of bins where 1) neither bin is in any region of interest AND 2) the distance between the bins is in the same range as all bin-pairs from the regions of interest.
For convenience specify titles + colors for all possible pairs of regions to compare for plotting
```{r colors_and_titles}
# Map colors and titles to each region of contacts
colors.and.titles <- 
    tribble(
    ~region, ~region.title, ~region.color,
    'chr16p.deletion vs chr16p.telomere', '16p11.2 CNV - Telomere',                     '#526ab1',
    'chr16p vs chr16p',                   '16p.Distance-matched control',               '#e32528',
    'chr16q vs chr16q',                   '16q.Distance-matched control',               '#5a0c10',
    'chr16p vs chr16q',                   '16p-q.Distance-matched control',             '#9e1517',
    'chr16p.deletion vs chr16p',          '16p11.2 CNV - 16p.Distance-matched control', '#8000ff',
    'chr16p.deletion vs chr16q',          '16p11.2 CNV - 16q.Distance-matched control', '#bf80ff',
    'chr16p.telomere vs chr16p',          'Telomere - 16p.Distance-matched control',    '#330066'
    ) %>% 
    mutate(region.title=factor(region.title))
```
Final maual check that regions have the correct number of bin-pairs.
```{r manual_check_regions, eval=FALSE}
plot.df %>% 
    filter(region.title %in% c('16p.Distance-matched control', '16p11.2 CNV - Telomere')) %>% 
    group_by(Sample.ID, IF.Threshold, region) %>%
    mutate(group=cur_group_id()) %>%
    ungroup() %>%
    filter(group == 2) %>% 
    count(range1, range2) %>% print(n=Inf)
```

## Weiner et al. 2022 Style Analysis{.tabset}

Here I just recreate the figure+results as simply as possible

## Direct Replication

```{r basic_boxplot, results='asis', message=FALSE, fig.height=5, fig.width=7} 
plot.df %>% 
    filter(
        region.title %in% c(
            '16p.Distance-matched control',
            '16p11.2 CNV - Telomere'
        )
    ) %>% 
    make_nested_plot_tabs(
        group_cols=c('IF.Threshold', 'Celltype', 'isMerged'),
        max_header_lvl=3,
        plot_fnc=plot_contacts_regions_boxplot,
        bin_col='region.title',
        count_col='IF',
        color_col='region.color',
        facet_row='ReplicateNum',
        facet_col='Genotype',
        tip.length=0,
        expansion=c(0.05, 0.01, 0.1, 0.01),
        plot.margin=margin(0.05, 0.05, 0.05, 0.6, "in")
        n_size=3,
        test_offset_y=4,
        # n_pos=0,
        # n_size=2.5,
        # test_offset_y=7,
    )
```

## Include more regions

We have more data so we can visualize some other relevant comparisons for this CNV region and use corrected p-values.
```{r extra_regions, results='asis', message=FALSE, fig.width=11, fig.height=11, eval=FALSE} 
plot.df %>% 
    filter(
        region.title %in% c(
            '16p11.2 CNV - Telomere',
            '16p.Distance-matched control',
            '16p-q.Distance-matched control',
            '16q.Distance-matched control',
            '16p11.2 CNV - 16p.Distance-matched control',
            '16p11.2 CNV - 16q.Distance-matched control',
            'Telomere - 16p.Distance-matched control'
        )
    ) %>% 
    plot_contacts_regions_boxplot(
        bin_col='region.title',
        count_col='IF',
        color_col='region.color',
        facet_row='ReplicateNum',
        facet_col='Genotype',
        p.adjust.method='BH',
        p.max=0.05,
        tip.length=0,
        test_offset_y=3
    )
```

# Now Repeat for Other RGDs

## 8p23 CNV - Telomere Contacts

