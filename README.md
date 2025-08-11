# 16p HiC Analysis

Analysis and code of HiC data for 16p, NIPBL and WAPL cells.

## Links to figures/notebooks

NIPBL+WAPL
- /data/talkowski/Samples/WAPL_NIPBL/HiC/notebooks/hicrep.html
- /data/talkowski/Samples/WAPL_NIPBL/HiC/notebooks/matrix.QC.html
- /data/talkowski/Samples/WAPL_NIPBL/HiC/notebooks/TADs.html
- /data/talkowski/Samples/WAPL_NIPBL/HiC/results.iN/plots/

16p
- /data/talkowski/Samples/16p_HiC/notebooks/hicrep.html
- /data/talkowski/Samples/16p_HiC/notebooks/matrix.QC.html
- /data/talkowski/Samples/16p_HiC/notebooks/replicate.robinson.results.html

## Repo+Results structure

Analysis and code of HiC data for 16p editied cell lines.
Notice that all files contain a `SampleID` specifying which library the results are.
Some files contain a pair of `SampleID`s since they are comparing 2 samples e.g. HiCRep.
All SampleIDs are follow the same format format `Project.CellType.Genotype.BioRepID.TechRepID` e.g. 16p.NSC.DEL.A3.TR1 

- `16p`: project this sample is a part of
- `NSC`: Celltype {NSC, iN}
- `DEL`: Genotype of the sample for the region/gene of interest {WT,DEL,DUP}
- `A3`: ID string specifiying which biological replicate the sample is
- `TR1`: ID string specifiyng which technical replicate the sample is


```
# Note ... indicates that these files exist for all samples/pairs, truncated for brevity
./
├── public.data/             # Public HiC datasets to compare against
├── distiller-nf/            # distiller nextflow installation
├── notebooks/               # results+figures
├── scripts/                 # scripts+notebook backends
├── reference.files/         # mostly conda envs, coordinates etc.
│   ├── cooltools.env.yml
│   ├── distiller.env.yml
│   ├── TADLib.env.yml
│   └── README.md
├── 16p.sample.metadata.tsv     # metadata for all HiC samples
├── fastq/                  # Raw reads for HiC samples
│   ├── 22LCC2LT4_3_2148261314_16pDELA3NSCHiC_S1_L003_R1_001.fastq.gz
│   ├── 22LCC2LT4_3_2148261314_16pDELA3NSCHiC_S1_L003_R2_001.fastq.gz
│   └── ....fastq.gz
├── sample.configs/         # Config files for distiller
│   ├── 16p.NSC.DEL.A3.TR1.distiller.yml
│   └── ....distiller.yml
└── results/                # all HiC results
    ├── sample.QC/
    │   ├── multiqc.reports/
    │   │   ├── fastp.multiqc.html
    │   │   ├── fastqc.multiqc.html
    │   │   ├── pairtools.multiqc.html
    │   │   └── qc3C.multiqc.html
    │   └── qc3C/
    │       ├── 16p.NSC.DEL.A3.TR1/    # SampleID, unique to each sample (library)
    │       │   └── report.qc3C.json 
    │       └── ../
    ├── coolers_library/
    │   ├── 16p.NSC.DEL.A3.TR1/
    │   │   ├── 16p.NSC.DEL.A3.TR1.hg38.mapq_30.1000.cool
    │   │   ├── 16p.NSC.DEL.A3.TR1.hg38.mapq_30.1000.mcool
    │   │   ├── 16p.NSC.DEL.A3.TR1.hg38.no_filter.1000.cool
    │   │   └── 16p.NSC.DEL.A3.TR1.hg38.no_filter.1000.mcool
    │   └── .../
    ├── fastqc/
    │   ├── 16p.NSC.DEL.A3.TR1/
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.0.2_fastqc.html
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.0.2_fastqc.zip
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.0.1_fastqc.html
    │   │   └── 16p.NSC.DEL.A3.TR1.lane1.0.1_fastqc.zip
    │   └── .../
    ├── mapped_parsed_sorted_chunks
    │   ├── 16p.NSC.DEL.A3.TR1/
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.hg38.0.fastp.json
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.hg38.0.fastp.html
    │   │   ├── 16p.NSC.DEL.A3.TR1.lane1.hg38.0.pairsam.gz
    │   │   └── 16p.NSC.DEL.A3.TR1.lane1.hg38.0.bam
    │   └── .../
    └── pairs_library
        ├── 16p.NSC.DEL.A3.TR1/
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.dedup.stats
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.dups.bam
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.dups.pairs.gz
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.nodups.bam
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.nodups.pairs.gz
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.unmapped.bam
        │   ├── 16p.NSC.DEL.A3.TR1.hg38.unmapped.pairs.gz
        │   └── 16p.NSC.DEL.A3.TR1.hg38.nodups.pairs.gz.px2
        └── .../
```

## Producing Results

### Running distiller pipeline

Each sample as a `.yml` file (in `./sample.configs`) specifiying the params for `distiller-nf` to run the sample with. We separte each sample into its own file so we can run them in parallel, but the only difference between files are the input fastq files, all pipeline parameters are the same. 

```bash
# launch 1 slurm job running the distiller-nf pipeline per sample 
$ ./scripts/run.distiller.sh             
        -a $HOME/miniforge3              # location of conda install
        -w ./work                        # nextflow working dir
        -p bigmem                        # slurm --partition
        -m 40                            # slurm --mem (in Gb)
        -t 8                             # slurm --ntasks-per-node 
        -c 8                             # slurm --cpus-per-task
        ./sample.configs/*.distiller.yml # distiller-nf param files per sample
# copy pastable
./scripts/run.distiller.sh -a $HOME/miniforge3 -w ./work -p bigmem -m 40 -t 8 -c 8 ./sample.configs/*.distiller.yml
```
This is 64 cores total, with 128 maxCPUs set in the nextflow config. 
This fully processes a sample (`.fastq -> .mcool`) with ~400M reads in ~9h.

### Running qc3C profiling

Use the tool `qc3C` [github](https://github.com/cerebis/qc3C) in bam mode to profile quality metrics for our HiC samples.

```bash
$ ./scripts/matrix.utils.sh qc3c                     # run qc3C with specified params
        ./results/sample.QC/qc3C/                    # output dir
        DpnII HinfI                                  # enzymes used
        results/mapped_parsed_sorted_chunks/**/*.bam # bam files produced by distiller for each sample
```

### Generate MultiQC reports

distiller-nf outputs multiple QC reports/files per sample than can each be aggregated into a single multiqc report [docs](https://docs.seqera.io/multiqc). 
In the two previous steps we have also generated stats files which can be aggretaed via multiqc into their own reports. 
Ultimately we can generate 3 multiqc reports 

- `fastqc`: Quality statistics for reads
- `fastp`: Trimming+MAPQ statistics for filtered reads that were aligned
- `qc3C`: Quality statistics from sub-sampled aligned reads
- `pairtools stats`: Summary statistics of processed HiC pairs

```bash
$ ./scripts/matrix.utils.sh multiqcs         # generate multiqc reports for different outputs from distiller,qc3C
        ./results/sample.QC/multiqc.reports/ # output dir
        ./results/                           # all multiqc data is under here in specific directories
```

### Calculate QC Metrics

Two things we want to check to QC matrix samples 
1. pair frequency by distance 
2. Cis/Trans pair frequency
3. minimum resolution as defined in [Rao et al. 2014](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4?cc=y%3D).

Metrics 1 and 2 are calcualted by `distiller-nf` and found in the `*.dedup.stats` files. 
For metric 3 and subsequent plots we need to calculate the per-bin coverage using [cooltools coverage ](https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-coverage).

```bash
$ ./scripts/matrix.utils.sh coverage
        ./results/sample.QC/coverage/
        ./results/coolers_library/**/*.mapq_30.1000.mcool
```

### Merging Matrices

For some subsequent steps we will analyze matrices formed by merging all biological/technical replicates for a given Edit+Genotype+CellType (e.g. 16p+WT+NSC).
Merging her means summing all the total number of contacts over all samples for each bin-bin pair (matrix entry) and is handled by [cooler merge](https://cooler.readthedocs.io/en/latest/cli.html#cooler-merge).
We separately create merged matrices for MAPQ filtered and unfiltered contacts for each 

```bash
# Wrapper command to merge specific matrices for each genotype
$ ./scripts/matrix.utils.sh merge_16p_matrices ./results/pairsn
$ ./scripts/matrix.utils.sh merge_16p_matrices ./results/coolers_library
```
which produces the following files
```
results 
└── coolers_library
    ├── 16p.NSC.WT.Merged.TR1/
    │   ├── 16p.NSC.WT.Merged.TR1.hg38.mapq_30.1000.cool
    │   ├── 16p.NSC.WT.Merged.TR1.hg38.mapq_30.1000.mcool
    │   ├── 16p.NSC.WT.Merged.TR1.hg38.no_filter.1000.cool
    │   └── 16p.NSC.WT.Merged.TR1.hg38.no_filter.1000.mcool
    ├── 16p.NSC.DEL.Merged.TR1/
    │   └── ...
    └── 16p.NSC.DUP.Merged.TR1/
        └── ...
```

### HiCRep Analysis

We use HiCRep to calculate the "reproducibility score" for all pairs of sample matrices, under several parameter combinations. The command below actually runs the HiCRep and produces 1 file per sample pair + parameter combination, each file contains scores for each chromosome separately (`chr{1..22,X,Y}`).

```bash
# Compare all pairs of matrices where all contacts have both mates MAPQ > 30
$ ./scripts/run.hicrep.sh 
        ./results/hicrep/
        $(find ./results/coolers_library -maxdepth 99 -type f -name "*.mapq_30.1000.mcool")
# Compare all pairs of matrices with no MAPQ filtering
$ ./scripts/run.hicrep.sh 
        ./results/hicrep/
        $(find ./results/coolers_library -maxdepth 99 -type f -name "*.no_filter.1000.mcool")
```

After this there is a `.Rmd` notebook that coallates these files into a single neat dataframe that is used for plotting.

### TAD Analysis

We produce TAD annotations from the individual and merged matrices. 
We use the individual matrix annotations to assess (

#### TAD Annotation

We generate the TAD annotations we use 2 different programs:

1. [HiTAD](https://xiaotaowang.github.io/TADLib/domaincaller.html)
2. [cooltools insulation](https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-insulation)

For `HiTAD` we only generate single-level TADs (not hierarchical) using default parameters. This produces a set of bin-pairs, each pair marking the start and end of the predicted TAD.

For `cooltools insulation` it only annotates whether a bin is a TAD boundary or not, it does not group a pair of boundaries to explicitly define the start/end of a specific TAD.
This changes the downstream analysis, but it is easy to adjust for the sake of making comparisons between 2 different TAD annotations.
For cooltools we also compute TAD annotations for multiple sets of hyper-parameters, just to asses how much these parameters affect the annotations.

Both of these tools also output the Diamond Insulation (DI) score calcualted per genomic bin.
This can be used to score regions to compare the relative insulation of regions. 
This likely only makes sense to calculate between boundaries or maybe just compare the distribution of all insulation scores across all bins within 2 regions.
No specific plans to do this analysis yet
The DI annotation data can also be compared directly to ATAC-Seq data to look for regions of high correlation between ATAC-Seq tracks and DI tracks.
We can also check if is ATAC peaks and and insulation boundaries are close/overlap (midpoint distance?)

#### TAD Comparison

For 2 different TAD annotation sets (start/end pairs) of the same region (e.g. chr16) we can calculate the similarity of the 2 sets by computing the Measure of Concordance as defined in [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1596-9#Sec21). 

For 2 different TAD boundary annotation sets (just is/is not a boundary) of the same region we can calculate similarity as follows:
1. Calculate the distance between all pairs of boundaries between the sets
2. For each boundary in set 1 pick the boundary in set 2 with the samllest distance
3. If there are any boundaries are set 2 that are the nearest to 2 different set 1 bounadaries, assign it to the closest set 1 boundary and assign the 2nd closest boundary in set 2 to the remaining set 1 boundary
4. Using all the boundary-pair distance, summarize in at least 1 of the following ways
   1. Test if the distribution is > 0 (KS-test vs Gaussian + observed variance) -> need to correct p-values since there are many tests (1 per pair of TAD annotation set)
   2. Calculate mean distance and just compare that between pairs

### Loop Analysis

#### Loop Annotation

#### Loop Comparison

### Compartment Analysis

#### Compartment Annotation

#### Compartment Comparison

### Differential Contact Analysis

