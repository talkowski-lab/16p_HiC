# Notes on HiC Analysis

## Paper List

### All Papers

#### []()

#### []()

#### []()

#### [CTCF shapes chromatin structure and gene expression in health and disease](https://pmc.ncbi.nlm.nih.gov/articles/PMC9442299/)

#### [Cohesin and CTCF differentially affect chromatin architecture and gene expression in human cells](https://www.pnas.org/doi/full/10.1073/pnas.1317788111)

#### [Spatial patterns of CTCF sites define the anatomy of TADs and their boundaries](https://link.springer.com/article/10.1186/s13059-020-02108-x#Abs1)

#### [Epigenetic and 3D genome reprogramming during the aging of human hippocampus](https://www.biorxiv.org/content/10.1101/2024.10.14.618338v1.full.pdf)

#### [Comparison of computational methods for Hi-C data analysis](https://www.nature.com/articles/nmeth.4325)  

#### [Revisiting Assessment of Computational Methods for Hi-C Data Analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531246/#sec2-ijms-24-13814)  

#### [A comprehensive review and benchmark of differential analysis tools for Hi-C data](http://clementine.wf/doc/pdf/jorge_etal_BB2025.pdf)  

#### [A comparison of topologically associating domain callers over mammals at high resolution](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04674-2#Sec18)  

#### [Comparative study on chromatin loop callers using Hi-C data reveals their effectiveness](https://www.biorxiv.org/content/biorxiv/early/2023/11/25/2023.11.24.567971.full.pdf)  

#### [A comprehensive benchmarking with interpretation and operational guidance for the hierarchy of topologically associating domains](https://www.nature.com/articles/s41467-024-48593-7?error=cookies_not_supported&code=69e5a830-5d27-4ae4-b469-757ead4ae67e)  

#### [Comparison of computational methods for the identification of topologically associating domains](https://link.springer.com/article/10.1186/s13059-018-1596-9#Abs1)  

#### [Unfolding neural diversity: how dynamic three-dimensional genome architecture regulates brain function and disease](https://www.nature.com/articles/s41380-025-03056-3?fromPaywallRec=false)  

#### [TAD Calling walkthrough](https://zhonglab.gitbook.io/3dgenome/chapter2-computational-analysis/3.2-higer-order-data-analysis/tad-calling-algorithms)  

#### [IDR2D identifies reproducible genomic interactions](https://academic.oup.com/nar/article/48/6/e31/5721212)  

#### [HiC-Bench Differential TAD Boundaries](https://github.com/NYU-BFX/hic-bench/wiki/Differential-boundary-insulation)

#### [Removing unwanted variation between samples in Hi-C experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC11074651/)  

#### [Systematic evaluation of chromosome conformation capture assays](https://www.nature.com/articles/s41592-021-01248-7)  

#### [On the assessment of statistical significance of three-dimensional colocalization of sets of genomic elements](https://academic.oup.com/nar/article/40/9/3849/1138988)  

#### [Local and global chromatin interactions are altered by large genomic deletions associated with human brain development](https://pmc.ncbi.nlm.nih.gov/articles/PMC6297223/#Abs1)

* “lymphoblastoid cell lines (LCLs) derived from patients with 22q11DS, that chromatin marks, chromatin domains, and long-range chromosome interactions are affected in several distinct ways by the large, common, and strongly disease-associated CNV on chromosome 22q11.2”  
* “11 human LCLs (5 patient cell lines with 22q11.2 deletion and 6 control cell lines without), with a total of 3.1 billion Hi-C contact reads of which 680 million read-pairs were of high quality and used for the downstream analyses”  
* “These phased SNVs were then used to create haplotype-specific Hi-C interaction maps. Haplotype-specific Hi-C analysis requires very deep sequencing coverage of Hi-C libraries since only paired-end reads that cover at least one informative (heterozygous) SNV can be used for that purpose. To increase the Hi-C coverage for chromosome 22q to these required very deep levels, we carried out custom-designed chromosome-wide targeted capture Hi-C.”  
* “We found that not all of the available normalization methods are robust for the use with Hi-C data coming from genomes with large CNVs. However, the hicpipe algorithm22 is quite suitable for this purpose”  
* “This decrease is consistent with the copy number of the 22q11.2 deletion region in the patient cell lines, as all of the 22q11DS cell lines are heterozygously deleted for this region. No such decrease of chromosomal contacts that involved an extended and contiguous chromosomal region was observed that did not involve the 22q11.2 region.”  
* “ We found the trans-contacts involving the 22q11.2 deletion in the patient cell lines and any other chromosome also decreased compared to control cell lines”  
* “The chromosome contacts between proximal and distal flanking regions of the 22q11.2 deletion increased on the chromosome 22q with the deletion when compared with the intact chromosome 22q within the same patient cell line”  
* “We did not observe differences in the A/B compartments between cell lines from patients and from controls and neither between the two homologous chromosomes 22q in any individual cell line, patient, or control. Our results indicated that A/B compartments of the homologous chromosome with the 22q11.2 deletion were not affected by the deletion.”  
* “Although there were variations in the calling of topological domains across individuals and between the homologous chromosomes, the direction indices were highly consistent between the two homologous chromosomes in the controls (Supplementary Fig. 6c-d, Supplementary Fig. 7\) as well as in the patients (Fig. 4, Supplementary Fig. 6a-b), suggesting no changes of topological domains on the homologous chromosome with the 22q11.2 deletion.”  
* “The second largest fold change (1.96) of cis-contacts involving region 21.5–22 Mbp was for contacts with region 50–51 Mbp, i.e., toward the very telomeric end of chromosome 22q” \-\> validated with 3D FISH”  
* “Inter-chromosomal contacts of GM06990 as determined by our own Hi-C data for this line (Supplementary Fig. 2a) show the same patterns of chromosomal interactions across the nucleus as in ref. 20; i.e., small chromosomes generally have more interactions with each other than larger chromosomes with each other and many more than chromosomes in the medium size range.”  
* On the genome-wide level, we found 272 trans-contacts with a Fisher’s exact test p value of \<0.0001 (Fig. 6a). Interestingly the majority of these chromosomal trans-contacts did not involve chromosome 22q as one of the interacting partners.” \-\> reconfirmed with permutations  
* Very basic DEG \+ Functional Enrichment analysis  
* “As a cautionary note on the technical level, we demonstrated here that for genomes with a large deletion CNV the appropriate normalization methods for Hi-C data have to be chosen with great care to avoid false findings. For instance, we would have reached the conclusion that the chromosomal contacts within the deletion regions are not decreased in cell lines with deletion compared with control cell lines if the hiclib software package23 had been used for normalization of Hi-C data, instead of the hicpipe package22.” \-\> should visualize the ICE normalized contacts in our data for the deletion regions  
* METHODS  
  * Used iterative mapping  
  * “read pairs whose sum of distances from mapped positions to the nearest restriction sites is larger than the length of the fragments in the Hi-C library were further removed by hicpipe”  
  * Call TADs on 40Kb bins  
  * All other analyses 500Kb  
  * Differential trans contact analysis with `t.test()`   
  * Permuted genotype labels

#### [Transcription shapes 3D chromatin organization by interactingwith loop extrusion](https://www.pnas.org/doi/pdf/10.1073/pnas.2210480120)

* “We observed new contacts between cohesin islands that appeared as Hi-C “dots” (island–island dots) that bridged distant genomic sites, consistent with the formation of cohesin-mediated chromatin loops”  
* “Wapl KO increases cohesin residence time,thus increasing the number of cohesin complexes on chromatin (55) and allowing time for accumulation in islands”   
* Figure S1 shows pile-up patterns around CTCF sites for WT, WAPL KO![][image1]  
* “In WT and CTCF and Wapl mutants, this revealed that individual genes are insulating, and active genes generate stronger insulation than inactive genes”  
* “Insulation is also weakened in Wapl KO and DKO cells (Fig. 1B), where increased residence time presumably allows loop-extruding cohesins to traverse the gene and bring regions upstream and downstream of the gene into contact.  
* ![][image2]  
* 

#### [Comparison of computational methods for Hi-C data analysis](https://www.nature.com/articles/nmeth.4325)

* “To compare the reproducibility of interactions called in different replicates, we calculated the similarity coefficient of Jaccard Index (Jaccard Index, JI) as a measure of the overlap between sets of interactions.For most data sets, the reproducibility among replicates of the same data set (intra-data set) was low at all resolutions (Fig. 2d and Supplementary Fig. 3a), yet it was significantly higher than random sets of interactions (P values ≤ 0.001; Supplementary Fig. 3b).”  
* “The intra-data set reproducibility remained similar when comparing replicates of the same cell line processed using different restriction enzymes (Supplementary Fig. 5). However, the inter-data set reproducibility—i.e., the concordance between interactions called in samples of the same cell line in different data sets (using different protocols and enzymes)—was much lower (median JI \< 4 × 10−4; Supplementary Fig. 6).”

#### [Revisiting Assessment of Computational Methods for Hi-C Data Analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531246/#sec2-ijms-24-13814)

#### [Analysis methods for studying the 3D architecture of the genome](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)

#### [Unfolding neural diversity: how dynamic three-dimensional genome architecture fregulates brain function and disease](https://www.nature.com/articles/s41380-025-03056-3?fromPaywallRec=false)

#### [Rules of engagement for condensins and cohesins guide mitotic chromosome formation](https://www.biorxiv.org/content/biorxiv/early/2024/04/30/2024.04.18.590027.full.pdf)

#### [Epigenetic and 3D genome reprogramming during the aging of human hippocampus](https://www.biorxiv.org/content/10.1101/2024.10.14.618338v1.full.pdf)

#### [Topological domains in mammalian genomes identified by analysis of chromatin interactions](https://pmc.ncbi.nlm.nih.gov/articles/PMC3356448/)

#### [Condensin-Driven Remodeling of X-Chromosome Topology during Dosage Compensation](https://pmc.ncbi.nlm.nih.gov/articles/PMC4498965/#S1) 

#### [Topologically associating domain boundaries that are stable across diverse cell types are evolutionarily constrained and enriched for heritability](https://pmc.ncbi.nlm.nih.gov/articles/PMC7895846/)

#### [Spatial patterns of CTCF sites define the anatomy of TADs and their boundaries](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02108-x)

#### [Integrated Analysis of Hi-C and RNA-Seq Reveals the Molecular Mechanism of Autopolyploid Growth Advantages in Pak Choi (Brassica rapa ssp. chinensis)](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.905202/full)

#### [HiC-ACT: improved detection of chromatin interactions from Hi-C data via aggregated Cauchy test](https://www.cell.com/ajhg/fulltext/S0002-9297(21)00009-4)

####  [Rules of engagement for condensins and cohesins guide mitotic chromosome formation](https://www.biorxiv.org/content/biorxiv/early/2024/04/30/2024.04.18.590027.full.pdf)  

####  [Epigenetic and 3D genome reprogramming during the aging of human hippocampus](https://www.biorxiv.org/content/10.1101/2024.10.14.618338v1.full.pdf)  

####  [Topological domains in mammalian genomes identified by analysis of chromatin interactions](https://pmc.ncbi.nlm.nih.gov/articles/PMC3356448/)  

####  [Condensin-Driven Remodeling of X-Chromosome Topology during Dosage Compensation](https://pmc.ncbi.nlm.nih.gov/articles/PMC4498965/#S1)   

####  [Topologically associating domain boundaries that are stable across diverse cell types are evolutionarily constrained and enriched for heritability](https://pmc.ncbi.nlm.nih.gov/articles/PMC7895846/)  

####  [Spatial patterns of CTCF sites define the anatomy of TADs and their boundaries](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02108-x)  
####   [CTCF shapes chromatin structure and gene expression in health and disease](https://pmc.ncbi.nlm.nih.gov/articles/PMC9442299/)

#### Unformatted Paper Links

https://academic.oup.com/bib/article/23/1/bbab509/6446983#327092372
https://www.cell.com/ajhg/fulltext/S0002-9297(21)00009-4#fig1
https://www.nature.com/articles/nmeth.4325
https://pmc.ncbi.nlm.nih.gov/articles/PMC10531246/#sec2-ijms-24-13814
https://github.com/mdozmorov/HiC_data?tab=readme-ov-file#cell-lines
https://www.sciencedirect.com/science/article/pii/S2211124720309104?via=ihub
https://pmc.ncbi.nlm.nih.gov/articles/PMC6297223/#Abs1
https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.905202/full
https://www.nature.com/articles/s41380-025-03352-y.pdf
https://www.biorxiv.org/content/10.1101/2023.11.24.567971v1.full.pdf
https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full
https://academic.oup.com/nar/article/51/20/11142/7301278#424819438
https://academic.oup.com/nar/article/40/9/3849/1138988
https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad152/7319106#425603573
https://academic.oup.com/gigascience/article/9/1/giz158/5695848
https://pmc.ncbi.nlm.nih.gov/articles/PMC11413648/
https://pmc.ncbi.nlm.nih.gov/articles/PMC6408958/#S11
https://pmc.ncbi.nlm.nih.gov/articles/PMC6028237/
https://www.nature.com/articles/s41467-024-52296-4
https://www.nature.com/articles/s41594-024-01431-2?fromPaywallRec=false
https://www.nature.com/articles/s41592-021-01248-7#Fig2
https://www.nature.com/articles/nature11279
https://www.nature.com/articles/s41586-022-04570-y?fromPaywallRec=false
https://www.nature.com/articles/s41588-025-02489-4
https://www.nature.com/articles/s41380-025-03352-y#Sec8
https://www.nature.com/articles/s41467-020-19452-y
https://www.nature.com/articles/nature14222
https://www.nature.com/articles/nature19847#accession-codes
https://www.nature.com/articles/s41594-022-00892-7?error=server_error#Sec11
https://www.nature.com/articles/s41467-020-19283-x#Fig2
https://www.nature.com/articles/s41467-023-39043-x#Sec10
https://www.nature.com/articles/s41576-023-00638-1#Fig2

## QC Notes

* “However, between replicates, TAD structures are shown to share only of their boundaries, suggesting that chromosome structure is not a static feature, but remains variable even in identical cell populations” \- [discussion sectio](https://academic.oup.com/nargab/article/2/1/lqz008/5576144#209716725)n  
* Aligning Hi-C reads may further require adjusting aligner settings. Some aligners assume DNA libraries contain only contiguous fragments, leading to ‘mate rescue’ where one read’s alignment is modified or even forced based on its pair’s alignment. This behavior is incompatible with Hi-C, which produces chimeric molecules with unrelated alignments on each side. To avoid erroneous results, disable mate rescue/pairing and align reads pairs independently. [In bwa mem, use the ‘-SP’ flags to achieve this.](https://pairtools.readthedocs.io/en/latest/parsing.html#aligner-settings)  
* [Pairtools pair type table](https://pairtools.readthedocs.io/en/latest/formats.html#pair-types)  
* [cooler balancing](https://github.com/open2c/cooler/discussions/419#discussioncomment-9570329) is genome-wide Iterative Correction (matrix balancing)   
* [Default fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file#quality-filter) quality filtering thresholds  
* The orientation convergence distance reported by ***pairtools stats*** can also be used to remove all Hi-C byproducts from binned output conservatively: Such filtering is, however, typically unnecessary as cooler and cooltools by default ignore the first two diagonals in all computations. This filter is sufficient to remove by-products of 4bp-cutter Hi-C and Micro-C at resolutions \> \= 1kb. ([pairtools paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164))  
* However, we find restriction-based filters unnecessary for more recently published Hi-C and do not include them in the standard *pairtools* pipeline. First, in our tests of recently published Hi-C datasets \[[1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.ref001)\], the statistical properties of pairs located far from and close to restriction sites proved nearly the same ([S2a Fig](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.s002)). Second, we found that unligated pieces of DNA can be removed by a simpler filter against short-distance pairs, which can be calibrated using strand-oriented scalings \[[42](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.ref042)\] ([S2b Fig](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.s002)). For downstream analyses in cooler \[[5](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.ref005)\], such by-products are removed by dropping pairs of bins with separations below a cutoff distance, which corresponds to removing a few central diagonals of a binned contact matrix. Finally, the annotation of restriction sites becomes less accurate and unproductive for libraries generated with frequent and/or flexible cutters (e.g., DpnII, MboI, and DdeI), cocktails thereof, and impossible in restriction-free protocols, such as Micro-C \[[10](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.ref010)\] and DNase Hi-C \[[53](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164#pcbi.1012164.ref053)\]. ([pairtools paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012164))

## Tools Considered

[Github](https://github.com/mdozmorov/HiC_tools) list of HiC tools

### Misc

#### HiC1Dmetrics 

[Docs](https://h1d.readthedocs.io/en/latest/)  
[Paper](https://h1d.readthedocs.io/en/latest/)

Calculate 1D metrics over the genome from a contact matrix i.e. compute 1 value per genomic bin

#### IDR2D 

[Paper](https://academic.oup.com/nar/article/48/6/e31/5721212)  
[R Package](https://bioconductor.org/packages/release/bioc/html/idr2d.html)  
Irreproducible discovery rate for HiC data, need to read paper

### HiC Reproducibility/Matrix Similarity

#### HiCRep

Paper  
Github  
Computes distrace-stratified correlation coefficient for each chr for a pair of HiC matrices.

#### HiCNoiseMeasuerer 

Paper  
[Github](https://github.com/JRowleyLab/HiCNoiseMeasurer)  
Quantifies how noisy a hic matrix is by calculating the rowswise mean auto-correlation. Since adjacent bins should be highly correlated, especially at lower resolutions, low acf means noisier data. Can be done at different levels of downsampling from the original matrix to get depth vs signal relationship.

#### HiC-Spector

[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870694)  
[Github](https://github.com/gersteinlab/HiC-spector)  
For a pair of hic matrices, compute the Lapacian for each and calculate the distance between the eigenvectors. More similar matrices (i.e. bio replicates vs diff cell lines) show more similarity by this method.

### Correction

#### HiCorr

[Github](https://github.com/JinLabBioinfo/HiCorr)  
Bias correction 

#### HiConfidence

[Paper](https://academic.oup.com/bib/article/24/2/bbad044/7033301)  
[Github](https://github.com/victorykobets/HiConfidence)  
Python tool for eliminating biases from the Hi-C data by downweighting chromatin contacts from low-quality (low-coverage) Hi-C replicates. For each replicate, calculate differential matrix and then the confidence for each pixel as the inverse pixel-wise difference divided by their mean, raised to the power of a tunable parameter k (Figure 2A). Used for correction for replicates' confidence in calculating TAD boundaries, intra-TAD densities, improves replicate reproducibility (stratum-adjusted correlation coefficient). Compared with multiHiCcompare, aids in differential analysis, compartment and TAD detection. Applied to D. melanogaster S2 cells, [GSE200078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200078), processed with distiller, pairtools, cooltools.

Normalization 

The R package [normGAM](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6331-8) compared methods for VC, SCN, ICE, KR2 and NLD normalizations, using implementations in matlab and R HiTC.

[This paper](https://www.tandfonline.com/doi/pdf/10.2144/btn-2019-0105) is a good review and comparison of different normalization methods and how they impact contacts \+ TADs

#### Binless normalization

[Github](https://github.com/3DGenomes/binless)

### TAD Calling

#### diffDomain

Python CLI tool  
[Paper](https://www.nature.com/articles/s41467-024-44782-6)  
[Github](https://github.com/Tian-Dechao/diffDomain?tab=readme-ov-file)  
[Docs](https://github.com/Tian-Dechao/diffDomain/wiki/1.1-Usage)

#### DiffTAD

Python CLI tool  
[Paper](https://www.biorxiv.org/content/10.1101/093625v1.full.pdf)  
[Github](https://bitbucket.org/rzaborowski/differential-analysis/src/master/)

#### TADCompare

R Package  
[Paper](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00158/full)  
[Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/TADCompare/inst/doc/Input_Data.html)

Provides a method to calculate a "boundary score" for each genomic bin, so it provides its own annotation of TADs. Given the boundary score a "differential boundary score" (a z-score according to authors) can be calculated for all bins for a pair of samples.  

#### hicDifferentialTAD

CLI Tool  
[Documentation](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDifferentialTAD.html)  
Paper

#### Localtadsim

CLI tool (Go)  
[Paper](https://academic.oup.com/bioinformatics/article/34/13/i475/5045717#405422926)  
[Github](https://github.com/Kingsford-Group/localtadsim)

**Update: Tool is broken and not update for 7 years**

A tool that compares how similar a set of TAD annotations is to another. Basically it treats a set of TAD annotations as a clustering on the genomic bins (each TAD is a cluster of bins) and compares the clusterings between two samples using the Variation of Information (distance metric). It can work with TAD annotations produced by basically any tool. It also provides a method to define "chromosomal intervals with statistically significantly similar TAD structure" between two samples using a permutation based approach.  
The R package [vssHiC](https://rpubs.com/nshokran/vssHiC) implements calling TADs with TopDOm, SpectralTAD and HiCseg.

#### Insulation Score

There are a few ways to do this but basically you compute an “insulation score” for every bin in the genome. A smaller insulation score usually means that you have fewer contacts crossing over a given bin and vice-versa. So bins with low insulation scores are more likely to be TAD boundaries. There are also several methods to define which specific bins are TAD boundaries. `Cooltools` provide functions for both computing an insulation score and defining TAD boundaries.

[cooltools insulation](https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html)

#### HiSeg

R Package  
[Paper](https://academic.oup.com/bioinformatics/article/30/17/i386/199711)  
[CRAN Archive](https://cran-archive.r-project.org/web/checks/2024/2024-03-24_check_results_HiCseg.html)

This package was removed from CRAN for compliance reasons but otherwise still seems fine to use.

#### PSYCHIC

Python CLI  
[Github](https://github.com/dhkron/PSYCHIC)  
[Paper](https://www.nature.com/articles/s41467-017-02386-3)

#### TADBit

Python CLI  
[Paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005665)  
[Github](https://github.com/3DGenomes/TADbit)  
[Documentation](https://3dgenomes.github.io/TADbit/)

#### TopDom

R package  
[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4838359/)  
[Github](https://github.com/HenrikBengtsson/TopDom)

#### CaTCH

R Package  
[Paper](https://genome.cshlp.org/content/27/3/479.long)  
[Github](https://github.com/zhanyinx/CaTCH)

Manually installable R package only? Weird and also minimal documentation

#### CHDF

Cpp executable  
[Paper](https://link.springer.com/article/10.1007/s40484-015-0047-9#Sec1)  
[Manual](https://static-content.springer.com/esm/art%3A10.1007%2Fs40484-015-0047-9/MediaObjects/40484_2015_47_MOESM2_ESM.cpp)  
[Executable](https://static-content.springer.com/esm/art%3A10.1007%2Fs40484-015-0047-9/MediaObjects/40484_2015_47_MOESM2_ESM.cpp)

Could work but dont know anything about Cpp, so would require more effort to run

#### HiTAD

CLI tool  
[Paper](https://academic.oup.com/nar/article/43/15/7237/2414371)  
[Documentation](https://xiaotaowang.github.io/TADLib/install.html)

Works\! No parameters except for Normalization method since it uses and adaptive window size method

#### HiCKey

CLI Tool \+ R package  
[CLI Github](https://github.com/YingruWuGit/HiCKey)  
[R Github](https://github.com/YingruWuGit/HiCKeyR)  
[Paper](https://www.biomedcentral.com/epdf/10.1186/s12859-021-04113-8?sharing_token=R3bA3nc3HNmQB5U0SmOtYm_BpE1tBhCbnbw3BuzI2RPMKplMNIEF9F6MFCiHuaH837SJaLhEfsLQewisIEJIazSrspW1CKDF3snTHtIshJsIDbENb910Vkl20IvDbr2uEbx0FcKvit0AOmuyrYmMIg4iEVSh9U0-_YKyMmgJnII%3D)

#### Amartus

CLI tool  
[Download Page](https://www.cs.cmu.edu/~ckingsf/software/armatus/)  
[Paper](https://link.springer.com/article/10.1186/1748-7188-9-14#Abs1)  
[Github](https://github.com/kingsfordgroup/armatus)

Out of date and basically never used by anyone except the authors.

### Compartment Analysis

#### dcHiC

[Github](https://github.com/ay-lab/dcHiC)  
Differential compartment analysis of HiC

### Loop Calling

#### IRD2D

[Paper](https://academic.oup.com/nar/article/48/6/e31/5721212)  
[Github](https://github.com/kkrismer/idr2d/)

Using loops called in individual technical replicates we can check their replicability across technical \+ biological replicates to find the most reproducible loops per condition (i.e. celltype \+ genotype)

#### Differential Contacts

[This review](https://academic.oup.com/bib/article/26/2/bbaf074/8051526#507324567) paper suggests that multiHiCComapare is well-calibrated and it allows us to use replicate information

## Matrix Comparison Strategies

### Contact Matrix Comparisons

**NOTE:** All the following analyses will be done per chromosome.  
These are ways we can quantify/compare individual contact matrices

Given our 6 single-sample WT matrices and our 3 aggregated WT matrices we can compute the following metrics per sample

1. Number of contacts over all bin pairs  
2. Total/Average contacts per bin (distribution)  
3. % empty bins/concordance of empty bins

Given these metrics/distributions we can compare them between all pairs of matrices, for distributions we can use MWU or KS tests to assess the differences.

The following stats are computed per pair of samples, 

1. Contact Ratio: distribution of individual contact ratios (matrix entries) over all bin pairs (e.g. corrected WAPL1\_matrix / corrected WAPL2\_matrix \~= 1 everywhere if the matrices are similar)  
2. [hicrep score:](https://pmc.ncbi.nlm.nih.gov/articles/PMC5668950/) basically a correlation coefficient for HiC matrices  
3. [HiC-spector:](https://github.com/gersteinlab/HiC-spector) also calculated a “reproducibility score” for a pair of HiC matrices

[This tool](https://www.bioconductor.org/packages/release/bioc/vignettes/diffHic/inst/doc/diffHicUsersGuide.pdf) provides a method (5.4.1) to estimate CNV effects on HiC contacts and remove them. If we have CNV calls for these samples then its worth visualizing/testing if regions with CNVS show increased contacts and which regions those are.

I think we can use `cooltools pileup` command to look at coverage around deletion regions to verify deletions as normally done with IGV

[This paper](https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC6130916&blobtype=pdf) defines the distal-to-local ratio i.e. ratio of contacts within a \+/-3Mb window around site *i* with contacts \> 3Mb away. 3Mb is chosen iio–.o-as most CTCF mediated interactions are \< 3Mb. This can be computed using the tool [h1d metrics](https://h1d.readthedocs.io/en/latest/onesample.html), can be could for identifying spots targeted by NIPBL+WAPL deletion.

### Comparing Pile-Ups

Should compare pile-ups of ICEd \+ O/E counts (cooltools expected)  
What is the control/calibration that this is informative?

* Compare to DEL flanking region pileups?  
* Compare against subset of random snippets

Comparing pile-ups of 

1. 2 sets of features within each contact matrix (i.e. CTCF sites vs random sites)  
2. Same feature set between 2 matrices (e.g. CTCF sites between WT and DEL)

How to compare pile-up matrices? 

1. Just calculate total distance pixel-wise between groups compared \-\> test if distributions of pixel-wise totals across snippet is \~ 0  
2. KS-test that distribution across all pixels is different   
   1. Do 1 ks test for each pixel (group 1 vs group2 samples)  
   2. Do 1 KS-test per pair of snippets

Pile-Ups to compare

* CTCF Sites  
* Large SCZ/Autims CNV flanking regions  
* RNA-Seq DEGs  
* Highly/lowly active genes  (based on TPM)  
* ATAC Seq Peaks

### Differential Contacts

[This review](https://academic.oup.com/bib/article/26/2/bbaf074/8051526#507324567) paper suggests that multiHiCComapare is well-calibrated and it allows us to use replicates as well. It also allows for a flexible design and we can include covariates as desired.
