# TODO

## cmds to run

- [ ] knit matrix.QC.Rmd
- [ ] knit matrix.coverage.Rmd
- [ ] knit hicrep.Rmd
- [ ] knit weiner.replication.Rmd
- [ ] knit TADs.Rmd
- [ ] knit TADCompare.Rmd
- [x] knit loops.Rmd
- [ ] knit loop.integration.Rmd
- [ ] knit loop.reproducibility.Rmd
- [ ] knit knit multiHiCCompare.Rmd
- [ ] knit HiC.Annotated.Heatmaps.Rmd

## HiCRep

- [x] generate HiCRep results
- [x] plot HiCRep results
  - [x] compare hyper-params
  - [x] compare genotypes + celltypes
  - [x] plot individual chromosomes

## Replicate Weiner et al. 2022 16p11.2 analysis

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

## TADCompare analysis

- [x] Generate TADCompare results using
  - [x] hiTAD inputs
  - [x] ConsensusTADs inputs
  - [x] SpectralTADs inputs
- [ ] plot detected differences

## Loop analysis

- [x] Call loops
- [x] Plot loops
  - [x] plot # of loops
  - [x] plot size of loops
  - [x] histograms/density of loop q-values
  - [ ] plot merged matrix depth vs total number of loops per condition
  - [x] plot volcano plot of loops
- [ ] quantify loop nesting
  - [ ] for every bin count # of loops it is inside
  - [ ] build tree of loop nesting structure
  - [ ] for every loop annotated is nesting level (i.e. tree depth)
- [x] IDR2D analysis
  - [x] compare how hyper-params affect loop reproducibility
  - [ ] plot loop reproducibility for chosen hyper-param set
- [ ] cCRE Integration
- [ ] DEG analysis
  - [ ] ???
- [ ] expression analysis
  - [x] map genes to loops and loops to reproducible/irreproducible
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

## phasing analysis

- [ ] get input files
  - [ ] .vcf for genomic background
  - [x] genome reference
- [ ] incorporated variants into refererence
- [ ] modify alignment params 
- [ ] adjust pairtools params
- [ ] generate phased HiC `.mcool` files
- [ ] merged phased replicates 


