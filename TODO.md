# TODO

## Read Papers

- [x] Original Rao et al 2014 paper [link](https://www.sciencedirect.com/science/article/pii/S0092867414014974?via%3Dihub#app1)
- [x] ENCODE guidelinkes [link](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/)
- [ ] Bimal's paper
- [ ] Rachita's CC paper [link](https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full)
- [ ] CHESS paper [link](https://www.nature.com/articles/s41588-020-00712-y)
- [ ] Serkan's iN paper [link](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-024-00538-6)
- [ ] Find QC stats/thresholds for calling loops
- [ ] Find QC stats/thresholds for calling compartments

## Lab Meeting 

- [ ] make pipeline diagram
- [ ] animate pairtools parsing slide
- [ ] review how pairtools parsing works
- [ ] make fastq+fastp stats table
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] distance decay plots
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] Pair statistics tables (pairtools multiqc)
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] line plots of normalized vs 
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] Plot contact heatmaps of specific regions at min resolution
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] QC Figures
  - [ ] remake bar plot figures for 16p
  - [ ] remake bar plot figures for NIPBL+WAPL
- [ ] Calculate minimum viable resolution table (Raot et al. 2014)
  - [ ] 16p samples
  - [ ] NIPBL+WAPL Samples
- [ ] HiCRep
  - [ ] re-read method details
  - [ ] generate results with 5e6 windowsize
  - [ ] plot hyper-parameter differences
  - [ ] NIPBL+WAPL
    - [ ] compare genotypes
    - [ ] compare individual vs merged
    - [ ] compare Edits
    - [ ] 1 example of per-chromosome scores
  - [ ] 16p
    - [ ] compare genotypes
    - [ ] compare individual vs merged
    - [ ] compare 6 original vs 3 new samples
    - [ ] 1 example of per-chromosome scores
- [ ] TADs
- [ ] multiHiCCompare
  - [ ] generate results for minimum resolutions
  - [ ] plot in regions of interest
  - [ ] calcualte genome-wide summary stats

## R01 Materials

- [ ] Write methods
  - [x] distiller pipeline
  - [ ] add statistics
    - [ ] fastq reads
  - [ ] how QC barplot stats are calculated
  - [ ] HiCRep method
  - [ ] TAD calling strategies
    - [ ] cooltools
    - [ ] hiTAD single-level
    - [ ] compare TAD similarity
- [x] QC barplot
  - [x] add cutoff lines based on ENCODE thresholds (15%, 35%)
- [x] HiCRep Heatmap
  - all samples vs all samples
  - get gw-hicrep by taking avg of all chr scores
  - unmerged matrices only
- [ ] contact heatmaps
  - [ ] raw contacts
  - [ ] balanced contacts (ICE)
  - [ ] log2 ratio between matrices
- [ ] TAD figures
  - [ ] number of TADs
  - [ ] TAD length
  - [ ] compare TAD similarity

## Code

- [x] MiSeq QC analysis
  - [x] 16p iNs
  - [x] 16p NSCs
- [x] calculate minimum viable resolution with Rao et al. 2014 definition 
- [x] create table of genomic regions to plot
- [x] plot read pair orientation for pairs
- [x] plot pair-category statistics with raw number (not %)
- [x] Plot HiCRep per chromosome (fix spacing)
- [ ] plot comparing raw vs balanced matrices
- [ ] create coverage line plot function
- [ ] create contact heatmap function
- [x] Finish TAD annotation results parser
- [ ] implement TAD similarity metrics
- [ ] finish re-write of multHiCCompare results generator

## 16p

- [x] run distiller to generate matrices
  - [9/9] NSCs
  - [3/3] iNs
- [ ] merge + balance matrices per Genotype+Celltype
  - [3/3] NSCs
  - [0/1] iNs
- [ ] calculate coverage for all matrices
  - [9/9] NSCs
  - [0/3] iNs
  - [3/4] merged samples
- [ ] calcualte HiCRep results 
  - [9/9] NSCs
  - [2/3] iNs
  - [3/4] merged samples
- [x] Pair QC notebook
- [x] HICRep notebook
- [x] recreate elise results in notebook
  - individual sample boxplots
  - merged sample boxplots
  - raw @ 100Kb
  - include p-values + no correction
- [ ] TAD analysis
  - [ ] cooltools
    - [9/9] NSCs
    - [2/3] iNs
    - [3/4] merged samples
  - [ ] hiTAD
    - [9/9] NSCs
    - [2/3] iNs
    - [3/4] merged samples
  - [ ] Plot summary statistics
    - [ ] TADs per chromosome
    - [ ] TAD sizes
    - [ ] TAD start position correlation/distribution ?
- [ ] MultiHiCCompare results
  - [ ] generate sparse matrix input
    - [0/9] NSCs
    - [0/3] iNs
    - [0/4] merged samples
  - [ ] generate results
    - [0/9] NSCs
    - [0/3] iNs
    - [0/4] merged samples
  - [ ] plot p-value distribution of contacts
    - [ ] pvalue histogram
    - [ ] Q-Q plot?
  - [ ] Compare stat distributions per result
    - [ ] log(10) p-value sum 
    - [ ] log(10) p-value median + var
    - [ ] number of significant contacts
      - [ ] at different thresholds
      - [ ] per chromosome

## NIPBL+WAPL

- [x] WAPL+NIPBL fix qc3C thing
- [x] merge + balance matrices per Genotype+Celltype
  - [5/5] iNs
- [x] calculate coverage for all matrices
  - [12/12] iNs
  - [5/5] merged samples
- [x] calcualte HiCRep results 
  - [12/12] iNs
  - [5/5] merged samples
- [x] Pair QC notebook 
- [x] HiCRep notebook
- [ ] TAD analysis
  - [x] generate TADs (10Kb, 25Kb, 50Kb, 100Kb)
    - [x] cooltools
      - [12/12] iNs
      - [5/5] merged samples
    - [x] hiTAD
      - [12/12] iNs
      - [5/5] merged samples
  - [ ] Plot summary statistics
    - [ ] TADs per chromosome
    - [ ] TAD sizes
    - [ ] TAD start position correlation/distribution ?
- [ ] MultiHiCCompare results
  - [ ] generate sparse matrix input
    - [0/12] iNs
    - [0/5] merged samples
  - [ ] generate results
    - [0/12] iNs
    - [0/5] merged samples
  - [ ] plot p-value distribution of contacts
    - [ ] pvalue histogram
    - [ ] Q-Q plot?
  - [ ] Compare stat distributions per result
    - [ ] log(10) p-value sum 
    - [ ] log(10) p-value median + var
    - [ ] number of significant contacts
      - [ ] at different thresholds
      - [ ] per chromosome

## Notes

### Plots of Merged Samples

#### Plot Types

- coverage line plot 
  - line plot across each region 
  - 1 line per sample
- heatmap of bin-bin contacts 
  - log10 raw contacts
  - balanced contacts

#### Plot list

##### 16p 

- WT
  - chr16
  - chr16p arm
  - chr16p deletion region +/- 1MB
  - chr16p deletion region +/- 100Kb
- DEL
  - chr16
  - chr16p arm
  - chr16p deletion region +/- 1MB
  - chr16p deletion region +/- 100Kb
- DUP
  - chr16
  - chr16p arm
  - chr16p deletion region +/- 1MB
  - chr16p deletion region +/- 100Kb

##### NIPBL+WAPL

- NIPBL WT
  - chr5
  - chr5p arm
  - chr5p13.2 arm 
  - NIPBL deletion region +/- 1Mb
  - NIPBL deletion region +/- 100Kb
- NIPBL DEL
  - chr5
  - chr5p arm
  - chr5p13.2 arm 
  - NIPBL deletion region +/- 1Mb
  - NIPBL deletion region +/- 100Kb
- WAPL WT
  - chr10
  - chr10q arm
  - chr10q23.3 arm
  - WAPL deletion region +/- 1Mb
  - WAPL deletion region +/- 100Kb
- WAPL DEL
  - chr10
  - chr10q arm
  - chr10q23.3 arm
  - WAPL deletion region +/- 1Mb
  - WAPL deletion region +/- 100Kb
- NIPBL+WAPL WT
  - chr10
  - chr10q arm
  - chr10q23.3 arm
  - WAPL deletion region +/- 1Mb
  - WAPL deletion region +/- 100Kb
  - chr5
  - chr5p arm
  - chr5p13.2 arm 
  - NIPBL deletion region +/- 1Mb
  - NIPBL deletion region +/- 100Kb

### TAD analysis

- finish generating annotations
- Plot summary statistics 
  - Individual + merged samples
  - TADs per chromosome
  - TAD sizes
  - TAD start position correlation?
- Compare annotations
  - compute MQC for a pair of TADs
  - compare individual matrices vs merged
  - compare merged vs merged
- Differential TAD calling
  - [TADCompare](https://dozmorovlab.github.io/TADCompare/)
  - [TADreg](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04614-0)
  - [DiffDomain](https://www.nature.com/articles/s41467-024-44782-6#Sec10)

### MultiHiCCompare

- individual + merged samples
- also compare individual vs merged samples
  - potential false positives?
- plot p-value distribution of contacts
  - pvalue histogram
  - Q-Q plot?
- Compare stat distributions per result
  - log(10) p-value 
    - sum 
    - median + var
  - number of significant contacts
    - different thresholds
    - per chromosome

### Loop Calling

- ENCODE says >= 2B unique paired-end reads to be "loop resolution"
