# TODO

## Tasks

- [ ] multiHiCCompare results with new samples
  - [ ] 16p Samples
  - [x] NIPBL+WAPL Samples
- [x] finish running HiC samples
- [x] write code to compare cooltools insulation TAD boundary annotations
- [x] write code to load cooltools and HiTAD annotations separtely
- [x] validate FC direction
- [ ] fix manhattan lines locations
- [ ] NIPBL,WAPL LFC vs LFC off diagonals -> find loations
- [ ] plot insulation scores correlations 
  - [ ] between biological replicates ?
  - [ ] between conditions per chr, 
  - [ ] comapre vs MHC hits
- [ ] compare significant bin-paris with insulation scores
- [ ] RBFOX1 location signals
- [ ] redo stuff with new 16p HiC samples
  - [x] create multiQC reports
  - [ ] merge matrices
  - [ ] compute coverage for minimum resolution
  - [ ] redo QC notebook
  - [x] run HiCRep on new matrix pairs
  - [x] redo HiCRep notebook
  - [ ] redo Weiner et al. 2022 notebook
  - [x] generate TAD annotations for new samples
  - [ ] redo TAD analysis notebook
- [ ] read loop callers review paper
- [ ] pick and install loop calling tools
  - [ ] mustache
  - [ ] ???
- [ ] write wrapper script to  call loops 

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

### multiHiCCompare 

- [x] volcanco plots
  - [x] 16p11.2 Samples
  - [x] NIPBL+WAPL

- [ ] manhattan plots
  - [ ] 16p11.2 Samples
    - [ ] chr16p11.2
    - [ ] chr16p
    - [ ] chr16
  - [x] NIPBL+WAPL
    - [ ] chr10q23.1 (WAPL)
    - [ ] chr10q (WAPL)
    - [ ] chr10 (WAPL)
    - [ ] chr5p13.2 (NIPBL)
    - [ ] chr5p (NIPBL)
    - [ ] chr5 (NIPBL)

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

## Read Papers

- [x] Original Rao et al 2014 paper [link](https://www.sciencedirect.com/science/article/pii/S0092867414014974?via%3Dihub#app1)
- [x] ENCODE guidelines [link](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/)
- [ ] Bimal's paper
- [ ] Rachita's CC paper [link](https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full)
- [ ] CHESS paper [link](https://www.nature.com/articles/s41588-020-00712-y)
- [ ] Serkan's iN paper [link](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-024-00538-6)
- [ ] Find QC stats/thresholds for calling loops
- [ ] Find QC stats/thresholds for calling compartments

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

