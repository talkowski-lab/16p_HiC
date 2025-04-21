# 16p HiC Analysis

Analysis and code of HiC data for 16p editied cell lines

File tree
```
# Note ... indicates that thee files exist for all samples, truncated for brevity
./
├── publicData          # Public HiC data 
├── distiller-nf        # distiller nextflow installation
├── notebooks           # results+figures
├── scripts             # scripts+notebook backends
├── reference.files     # mostly conda envs, coordinates etc.
│   ├── cooltools.env.yml
│   ├── distiller.env.yml
│   ├── TADLib.env.yml
│   └── README.md
├── sample.metadata.tsv # metadata for all HiC samples
├── fastq               # Raw reads for HiC samples
│   ├── 22LCC2LT4_3_2148261314_16pDELA3NSCHiC_S1_L003_R1_001.fastq.gz
│   ├── 22LCC2LT4_3_2148261314_16pDELA3NSCHiC_S1_L003_R2_001.fastq.gz
│   └── ....fastq.gz
├── sample.configs      # Config files for distiller
│   ├── 16p.DELA3.NSC.HiC.distiller.yml
│   └── ....distiller.yml
└── results.NSC         # all HiC results
    ├── sample.QC
    │   └── multiqc.reports 
    │       ├── fastp.multiqc.html
    │       ├── fastqc.multiqc.html
    │       ├── pairtools.multiqc.html
    │       └── qc3C.multiqc.html
    ├── coolers_library
    │   ├── 16p.DEL.A3.NSC.HiC/
    │   │   ├── 16p.DEL.A3.NSC.HiC.hg38.mapq_30.1000.cool
    │   │   ├── 16p.DEL.A3.NSC.HiC.hg38.mapq_30.1000.mcool
    │   │   ├── 16p.DEL.A3.NSC.HiC.hg38.no_filter.1000.cool
    │   │   └── 16p.DEL.A3.NSC.HiC.hg38.no_filter.1000.mcool
    │   └── .../
    ├── fastqc
    │   ├── 16p.DEL.A3.NSC.HiC/
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.0.2_fastqc.html
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.0.2_fastqc.zip
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.0.1_fastqc.html
    │   │   └── 16p.DEL.A3.NSC.HiC.lane1.0.1_fastqc.zip
    │   └── .../
    ├── mapped_parsed_sorted_chunks
    │   ├── 16p.DEL.A3.NSC.HiC/
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.hg38.0.fastp.json
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.hg38.0.fastp.html
    │   │   ├── 16p.DEL.A3.NSC.HiC.lane1.hg38.0.pairsam.gz
    │   │   └── 16p.DEL.A3.NSC.HiC.lane1.hg38.0.bam
    │   └── .../
    └── pairs_library
        ├── 16p.DEL.A3.NSC.HiC/
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.dedup.stats
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.dups.bam
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.dups.pairs.gz
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.nodups.bam
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.nodups.pairs.gz
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.unmapped.bam
        │   ├── 16p.DEL.A3.NSC.HiC.hg38.unmapped.pairs.gz
        │   └── 16p.DEL.A3.NSC.HiC.hg38.nodups.pairs.gz.px2
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
### Merging Matrices

For some subsequent steps we will analyze matrices formed by merging all biological/technical replicates for a given Edit+Genotype+CellType (e.g. 16p+WT+NSC).
Merging her means summing all the total number of contacts over all samples for each bin-bin pair (matrix entry) and is handled by [cooler](https://cooler.readthedocs.io/en/latest/cli.html#cooler-merge).
We separately create merged matrices for MAPQ filtered and unfiltered contacts for each 

```bash
# Wrapper command to merge specific matrices for each genotype
$ ./scripts/matrix.utils.sh merge_16p 
```
which produces the following files
```
results.NSC 
└── coolers_library
    ├── 16p.WT.Merged.NSC.HiC/
    │   ├── 16p.WT.Merged.NSC.HiC.hg38.mapq_30.1000.cool
    │   ├── 16p.WT.Merged.NSC.HiC.hg38.mapq_30.1000.mcool
    │   ├── 16p.WT.Merged.NSC.HiC.hg38.no_filter.1000.cool
    │   └── 16p.WT.Merged.NSC.HiC.hg38.no_filter.1000.mcool
    ├── 16p.DEL.Merged.NSC.HiC/
    │   ├── 16p.DEL.Merged.NSC.HiC.hg38.mapq_30.1000.cool
    │   ├── 16p.DEL.Merged.NSC.HiC.hg38.mapq_30.1000.mcool
    │   ├── 16p.DEL.Merged.NSC.HiC.hg38.no_filter.1000.cool
    │   └── 16p.DEL.Merged.NSC.HiC.hg38.no_filter.1000.mcool
    └── 16p.DUP.Merged.NSC.HiC/
        ├── 16p.DUP.Merged.NSC.HiC.hg38.mapq_30.1000.cool
        ├── 16p.DUP.Merged.NSC.HiC.hg38.mapq_30.1000.mcool
        ├── 16p.DUP.Merged.NSC.HiC.hg38.no_filter.1000.cool
        └── 16p.DUP.Merged.NSC.HiC.hg38.no_filter.1000.mcool
```
### Running qc3C profiling

Use the tool `qc3C` [github](https://github.com/cerebis/qc3C) in bam mode to profile quality metrics for our HiC samples.

```bash
$ ./scripts/matrix.utils.sh qc3c                             # run qc3C with specified params
        ./results.NSC/sample.QC/qc3C/                        # output dir
        results.NSC/mapped_parsed_sorted_chunks/16p.**/*.bam # bam files produced by distiller for each sample
# copy pastable
./scripts/matrix.utils.sh qc3c ./results.NSC/sample.QC/qc3C/ results.NSC/mapped_parsed_sorted_chunks/16p.**/*.bam
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
$ ./scripts/matrix.utils.sh multiqcs             # generate multiqc reports for different outputs from distiller,qc3C
        ./results.NSC/sample.QC/multiqc.reports/ # output dir
        ./results.NSC/                           # all multiqc data is under here in specific directories
# copy pastable
./scripts/matrix.utils.sh multiqcs ./results.NSC/sample.QC/multiqc.reports/ ./results.NSC/
```
### Restriction Fragment Analysis

We use `pairtools restrict` [docs](https://pairtools.readthedocs.io/en/latest/examples/pairtools_restrict_walkthrough.html) to annotate restriction fragments for each HiC read and use this for quality control.

The fist step is to generate the digested reference bed file for our data. This script Assumes genome fasta + bwa index are in the directory `${REF_DIR}`, specified as a variable in `./scripts/matrix.utils.sh`. 

```bash
# Generate multi-digested reference since we use ARIMA kit for HiC 
$ ./scripts/matrix.utils.sh digest_genome
# Annotate restriction fragments to reads with pairtools
$ ./scripts/matrix.utils.sh restrict                        # run pairtools restrict
        ./results.NSC/sample.QC/restriction.analysis/       # output dir
        ./results.NSC/pairs_library/16p.*/*.nodups.pairs.gz # pairsfiles for each sample
# copy pastable
./scripts/matrix.utils.sh restrict ./results.NSC/sample.QC/restriction.analysis/ ./results.NSC/pairs_library/16p.*/*.nodups.pairs.gz
```

### HiCRep Analysis

We use HiCRep to calculate the "reproducibility score" for all pairs of sample matrices, under several parameter combinations. The command below actually runs the HiCRep and produces 1 file per sample pair + parameter combination, each file contains scores for each chromosome separately (`chr{1..22,X,Y}`).

```bash
# Compare all pairs of matrices where all contacts have both mates MAPQ > 30
$ ./scripts/run.hicrep.sh 
        ./results//hicrep
        $(find ./results/coolers_library -maxdepth 99 -type f -name "*.mapq_30.1000.mcool" | grep -v "Merged")
# Compare all pairs of matrices with no MAPQ filtering
$ ./scripts/run.hicrep.sh 
        ./results/hicrep
        $(find ./results/coolers_library -maxdepth 99 -type f -name "*.no_filter.1000.mcool" | grep -v "Merged")
```

After this there is a `.Rnw` notebook that coallates these files into a single neat dataframe that is used for plotting.
