<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [TODO](#todo)
  - [Cohesin Results](#cohesin-results)
  - [cmds to run](#cmds-to-run)
  - [HiCRep](#hicrep)
  - [Replicate Weiner et al. 2022 16p11.2 analysis](#replicate-weiner-et-al-2022-16p112-analysis)
  - [TAD analysis](#tad-analysis)
  - [TADCompare analysis](#tadcompare-analysis)
  - [Loop analysis](#loop-analysis)
  - [multiHiCCompare analysis](#multihiccompare-analysis)
  - [Compartment analysis](#compartment-analysis)
  - [16p gghic Plots](#16p-gghic-plots)
  - [H1D analysis](#h1d-analysis)
  - [Gene Desert analysis](#gene-desert-analysis)
  - [TRADE analysis](#trade-analysis)
  - [Phasing analysis](#phasing-analysis)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# TODO

## Notebooks to Compile

- [x] knit Matrix.QC.Rmd
- [x] knit Matrix.Coverage.Rmd
- [x] knit HiCRep.Rmd
- [ ] knit Weiner.Replication.Rmd
- [ ] knit TADs.Rmd
- [ ] knit TADCompare.Rmd
- [x] knit Loops.Rmd
- [x] knit Loop.Reproducibility.Rmd
- [ ] knit Loop.Integration.Rmd
- [ ] knit multiHiCCompare.Rmd
- [ ] knit Compartments.Rmd
- [ ] knit HiC.Annotated.Heatmaps.Rmd

## HiCRep

- [x] generate HiCRep results
- [x] plot HiCRep results
  - [x] compare hyper-params
  - [x] compare genotypes + celltypes
  - [x] plot individual chromosomes
  - [x] similarity heatmaps
  - [x] compare tech reps to bio reps

## Replicate 16p11.2 analysis from Weiner et al. 2022

- [ ] wrap figures into notebook
- [x] do on merged matrices
  - [x] NSCs
  - [x] iNs
  - [x] compare genotypes for each region
- [-] do on individual matrices?
  - [-] NSCs
  - [-] iNs

## TAD analysis

- [x] Generate TAD Calls
  - [x] hiTAD
  - [x] ConsensusTAD using TADCompare
- [x] plot # of TADs
- [x] plot size of TADs
- [ ] plot TAD MoCs 
  - [x] between conditions
  - [ ] between our data and public data
- [ ] plot Intra/Inter TAD contact ratios
- [ ] DEG Integration
  - [ ] are DEGs in strong TADs?
  - [ ] are DEGs in more disimmilar TADs (MoC)
  - [ ] are DEGs in more disimmilar TADs (TADCompare)
  - [ ] are DEGs closer to boundaries than midpoints?
  - [ ] do DEGs often overlap TAD boundaries
- [ ] cCRE integration
  - [ ] is it more common for both anchors to be inside a TAD vs in different TADs?
  - [ ] compare size distributions as control

## TADCompare analysis

- [x] Generate TADCompare results using
  - [x] hiTAD inputs
  - [x] ConsensusTADs inputs
  - [x] SpectralTADs inputs
- [ ] plot detected differences
  - [ ] relative abundance of difference types
  - [ ] relative abundance of difference vs non-difference 

## Loop analysis

- [x] Call loops
- [x] Plot loops
  - [x] plot # of loops
  - [x] plot size of loops
  - [x] histograms/density of loop q-values
    - [ ] compare across conditions
      - [ ] between Celltype per Genotype
      - [ ] between Genotypes per Celltype
      - [ ] genome-wide
      - [ ] per chr
  - [ ] plot merged matrix depth vs total number of loops per condition
  - [x] plot volcano plot of loops
- [ ] APA analysis + comparison
  - [ ] 10% most significant loops
  - [ ] NR vs R loops between conditions
- [ ] check loop anchor association with CTCF sites
  - [ ] Permutation control?
  - [ ] other type of site comparison?
- [ ] quantify loop nesting
  - [ ] for every bin count # of loops it is inside
  - [ ] build tree of loop nesting structure
  - [ ] for every loop annotated is nesting level (i.e. tree depth)
- [x] IDR2D analysis
  - [x] compare how hyper-params affect loop reproducibility
  - [x] plot loop reproducibility for chosen hyper-param set
  - [ ] compare enrichment/qvalue of reproduced vs non-reproduced loops
- [ ] cCRE Integration
- [ ] DEG analysis
  - [ ] boxplot of logFC by condition (color by log pvalue)
- [ ] expression analysis
  - [ ] boxplot of expression by condition

## multiHiCCompare analysis

- [x] Generate multiHiCCompare results
- [x] plot DACs
  - [x] plot # of differential contacts 
  - [x] GW boxplots
  - [x] volcano plots
  - [x] manhattan plots
  - [x] log(FC) vs log(FC) plots, colored by log(GW adjusted p-value)

## Compartment analysis

- [ ] generate compartment calls
- [ ] try different bining strartegies
  - [ ] binary
  - [ ] quantile
- [ ] saddle plot of compartment interactions

## 16p gghic Plots

- [ ] wrap figures into notebook
- [x] regions
  - [x] chr16p11.2
  - [x] chr16 telomere
  - [x] chr1611.2 CNV 
  - [x] chr16p
  - [-] chr16q
  - [-] chr22q
  - [-] chr22q del region?
- [x] group samples
  - [x] plot all Genotypes per Celltype
  - [ ] plot all Celltypes per Genotype
- [ ] color loops by significance
- [ ] add gene annotations
- [ ] add track annotations 
  - [ ] loop nesting lvl (i.e. how many loops overlap each bin)
  - [ ] compartment phased eigenscore
  - [ ] marginal contact density
  - [ ] insulation score
  - [ ] gene desnity

## H1D analysis

- [ ] generate H1D metrics 

## Gene Desert analysis

- [ ] compare desert vs non-desert regions
  - [ ] loop density or total enrichment
  - [ ] TAD Density or total TAD 
  - [ ] mean/sum of TAD intra/inter contact ratios

## TRADE analysis

- [x] generate TRADE results from DESeq2 results
- [ ] Get shrunken logFCs from TRADE model 
- [ ] plot shrunken vs original logFCs 
- [ ] compare p.adj of TRADE DEGs vs rest

## Phasing analysis

- [ ] get input files
  - [ ] .vcf for genomic background
  - [x] genome reference
- [ ] incorporated variants into refererence
- [ ] modify alignment params 
- [ ] adjust pairtools params
- [ ] generate phased HiC `.mcool` files
- [ ] merged phased replicates 

