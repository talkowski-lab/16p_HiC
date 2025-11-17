# TODO

- list of sample group comparisons
  - iNs
    1. iN  DEL vs iN  WT
    1. iN  DUP vs iN  WT
    1. iN  DUP vs iN  DEL
  - NSCs
    1. NSC DEL vs NSC WT
    1. NSC DUP vs NSC WT
    1. NSC DUP vs NSC DEL
  - cross celltypes
    1. NSC WT  vs iN  WT
    1. NSC DEL vs iN  DEL
    1. NSC DUP vs iN  DUP

## Tasks

- [x] run qc3C on new samples
- [x] generate multiQCs on new samples
- [x] merge matrices
- [x] calculate bin-wise coverage on new samples
- [-] calculate expected coverage on new samples
- [ ] calculate HiCRep corr with new samples
  - [x] individual vs individual
  - [ ] merged vs merged
- [ ] generate TAD Results with merged matrices 
  - [ ] TAD calls with hiTAD
  - [ ] TAD boundary calls with cooltools
  - [ ] TADCompare results with spectralTAD calls
- [ ] generate multiHiCCompare results 
- [ ] generate Compartment results with merged matrices
  - [ ] dcHiC
  - [ ] HiCDOC
- [ ] generate H1DMetrics
  - individual and merged matrices
  - [ ] metrics
    - [ ] 
  - inferences
    - [ ] dTADs
    - [ ] stripeTAD
    - [ ] Hubs
- [ ] generate loop results
  - [ ] tools
    - [ ] cooltools dots
    - [ ] cLoops
    - [ ] mustache
  - [ ] generate pile-ups from loops

## Plots to Make

### Sample QC

- [x] barplot Pair Category bar plot per Sample
- [x] lineplot Pair Orientation and Frequency by Distance per Sample
- [x] heatmap Contacts between all chromosomes per Sample
- [x] contacts summary stats per samples per chr 
- [ ] lineplot plot comparing bin-wise coverage of balanced vs raw IFs

### HiCRep

- [x] boxplot of hyper-parameter comparison 
- [x] boxplot of all sample pairs scores per celltype per genotype
- [x] boxplot of all sample pairs scores per chromosome

### TADs

- [ ] hiTAD 
  - [x] boxplot number of TADs per Genotype per Chromosome
  - [x] boxplot TAD lengths per Genotype per Chromosome
  - [x] boxplot of MoCs between Genotypes
  - [x] heatmap mean MoC across pairs per Genotype per Chromosome

- [ ] cooltools 
  - [ ] boxplot number of boundaries per Genotype per Chromosome
  - [ ] boxplot number of boundaries per Parameter Set
  - [ ] boxplot boundary strengths per Genotype per Chromosome
  - [ ] boxplot boundary strengths per Parameter Set
  - [ ] boxplot inter-boundary lengths per Genotype per Chromosome
  - [ ] boxplot boundary strengths per Parameter Set
  - [ ] boxplot pairwise boundary concordance??? per Genotype per Chromosome
  - [ ] boxplot pairwise boundary strengths per Parameter Set

- [ ] TADCompare
  - [ ] ???

### multiHiCCompare 

- [x] boxplot of LFC + pvalue
- [x] volcanco plots
- [x] manhattan plots
- [x] plot fold-changes against eachother across comparisons
- [x] common DACs across conditions

### Loops

- [ ] cooltools dots results
  - [ ] number of loops per Param Set
  - [ ] number of loops per Genotype per Chromosome
  - [ ] APA plots of loops
  - [ ] distribution of loops/significance across each chr

- [ ] mustache loop results
  - [ ] number of loops per Param Set
  - [ ] number of loops per Genotype per Chromosome
  - [ ] APA plots of loops
  - [ ] distribution of loops/significance across each chr

