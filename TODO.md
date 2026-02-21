# TODO

- [x] validate FC direction
- [ ] fix manhattan lines locations

- [x] properly de-dup TAD MoC pairs used for correlation
- [ ] plot insulation lineplots
- [ ] plot insulation scores correlations

- [ ] NIPBL,WAPL LFC vs LFC off diagonals -> find locations
- [ ] compare significant bin-paris with insulation scores
- [ ] RBFOX1 location signals

- [?] generate GC track phase files for compartment calling
- [ ] look at other compartment calling tools
- [ ] generate compartment annotations

- [x] scb merge datasets
- [x] re-write to specify marker gene list per annotation+dataset
  - [x] make all the comparison directories
  - [x] write the data-set specific marker files there

- [x] rename RAD21 results with new genotypes
- [x] RAD21 MiSeq
- [x] CTCF MiSeq
- [x] make multiQC reports
- [x] produce coverage files

# Matrix QC + coverage

- [x] make Merged Matrices
  - WAPL.iN.DEL
  - WAPL.iN.WT
  - NIPBL.iN.DEL
  - NIPBL.iN.WT
  - All.iN.DEL
  - All.iN.WT

- [x] balance Merged Matrices

- [ ] calculate matrix bin-wise coverage
  - [x] individual
  - [ ] merged

- [ ] print minimum usable resolution per matrix 
  - [x] individual 
  - [ ] merged

- [x] heatmaps showing samples vs chr of 
  - [x] total contacts
  - [x] median contacts
  - [x] number of bin-pairs > 0
  - [x] plot line plot bin-wise total contacts across 

- [-] plot lineplot number of bin-pairs > 0 vs chr length in bins
  - one line per sample

# Call TADs method

- [x] Dont bother with individual TAD calling, only call on merged matrices

- Plot TAD comparisons
  - [x] heatmap N TADs         per condition per chr
  - [x] boxplot TAD Lengths    per condition
  - [x] boxplot TAD Lengths    per condition per chr
  - [x] boxplot TAD MoCs       per condition 
  - [x] boxplot mean TAD MoCs  per condition per chr
  - [x] heatmap mean TAD MoCs  per condition per chr

- per TAD metrics to compare
  - TAD intensities?
    - [ ] intra.TAD.contacts = sum of all contacts for all bin-pairs with BOTH   bins inside the TAD
    - [ ] outra.TAD.contacts = sum of all contacts for all bin-pairs with ONLY 1 bin  inside the TAD
    - normalize by
      - [ ] total number of bin the TAD covers (linear)
      - [ ] total number of bin-pairs with > 0 IF within the TAD
    - can compare distributions or summary stats
      - mean, median, min, max, 
      - Can do a t.test of intra vs intra between conditions
      - [ ] boxplot distribution of each TAD-pair differences
        - intra.TAD.contacts difference per TAD pair

- compare individual pairs of TADs?
  - Define the pairs of TADs to analyze i.e. the "same" TAD in both conditions
    - Not all TADs will have a pair
      - count how many TADs are unpaired after the matching
      - nTADs.unpaired.P1, nTADs.unpaired.P2
    - 2 TADs are a pair if they have the max "similarity" of all other TADs
    - Can define this for all TADs i.e. build a similarity matrix
      - cacluate all individual TAD similarities (MoC modified Jacard) 
      - filter all entries with similarity == 0
      - Solve the "Assignment Problem" using filtered similarity table (i.e. an edge-list)
      - The matching is all your pairs, ignore unpaired TADs

  - compare coordinates? 
    - total difference between start/end positions -> low total distance == same positions
      - [ ]    start.dist      = start.Numerator      - start.Denominator
      - [ ]      end.dist      =   end.Numerator      -   end.Denominator
      - [ ]     union.coverage = max(end.P2, end.P11) - min(start.P1, start.P2)
      - [ ] intersect.coverage = min(end.P2, end.P11) - max(start.P1, start.P2)
      - [ ] position.dist.bp   =     (abs(start.diff) + abs(end.diff)) 
      - [ ] position.dist.bins =     position.diff.bp / resolution
      - [ ] MoC

  - compare per TAD metrics

- examine cooltools insulation results?
  - [ ] compare param.combos
    - [x] heatmap N Boundaries
    - [x] boxplot Boundary strength
  - [ ] compare conditions
    - [x] heatmap N Boundaries      per condition per chr
    - [x] boxplot Boundary strength per condition 
    - [ ] boxplot Boundary strength per condition per chr
  - [ ] compare individual boundaries
    - position.dist = position.Numerator - position.Denominator
      - [ ] boxplot Min Boundary distances       per condition
      - [ ] boxplot Min Boundary distances       per condition per chr
    - strength.diff = strength.NUmerator - strength.Denominator
      - [ ] boxplot Min Boundary intensity diff  per condition
      - [ ] boxplot Min Boundary intensity diff  per condition per chr

- compare insulation profiles?
  - compare hyper.parameter differences
  - plot insulation across the genome
    - on chr5,10,1,X,17
    - one line per condition
  - calculate difference in insulation score
    - insluation.P1 - inslation.P2 
  - heatmap correlation Condition vs Condition
    - total insulation difference
      - per chr and genome-wide
      - boxplot of chrs per condition
MHC comparisons

# multiHiCCompare

- no merged comparison
- [x] individual matrices
  - [x] NIPBL.WT  vs  WAPL.WT  # background/noise/stochastic ? patient differences?
  - [x] NIPBL.DEL vs   All.WT  # NIBPL effect Main Result
  - [x] NIPBL.DEL vs NIBPL.WT  # NIBPL effect -> should replicate Main Result
  - [x]  WAPL.DEL vs   All.WT  #  WAPL effect Main Result
  - [x]  WAPL.DEL vs  WAPL.WT  #  WAPL effect -> should replicate Main Result
  - [x] NIPBL.DEL vs  WAPL.WT  # NIBPL effect -> should replicate Main Result
  - [x]  WAPL.DEL vs NIBPL.WT  #  WAPL effect -> should replicate Main Result
  - [-]   All.DEL vs   All.WT  
    - # ??? Common off targets? 
    - general deletion effect i.e. noise? -> no, it will just be all effects?
    - need to separate noise vs effects -> get intersection of NIPBL + WAPL deletion effects?
- use exact test
- [x] Plot significant contacts
- [x] Plot common hits across multiple comparisons
  - [x] NIPBL.WT vs WAPL.WT hits at the same significance lvl
  - [x] hits found in all 3 comparisons (robust to downsampling)
    - NIPBL.DEL vs   All.WT 
    - NIPBL.DEL vs NIBPL.WT
    - NIPBL.DEL vs  WAPL.WT

# TODO

## Tasks

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

### TADs

- [ ] plot insulation scores correlations 
  - [ ] between biological replicates ?
  - [ ] between conditions per chr, 

### Loops

- [ ] read loop callers review paper
- [ ] write wrapper script for loop calling
  - [x] cooltools dots
  - [ ] mustache
  - [ ] Fit-HiC-2

### multiHiCCompare

- [ ] multiHiCCompare results with new samples
  - [ ] 16p iN Samples
  - [x] Cohesin Samples
- [ ] fix manhattan lines locations
- [ ] Cohesin LFC vs LFC off diagonals -> find loations

### Integrations

- [ ] compare significant bin-paris with insulation scores
- [ ] RBFOX1 location signals

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

- [x] cooltools dots results
  - [x] number of loops per Genotype per Chromosome
  - [x] size of loops per Genotype per Chromosome
  - [x] distribution of loops/significance across each chr
  - [ ] APA plots of loops
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

