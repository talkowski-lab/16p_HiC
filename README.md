# 16p HiC Analysis

Analysis and code of HiC data for 16p editied cell lines

File tree
```
./
├── distiller-nf                    # installed distiller nextflow pipeline
├── fastq                           # fastq files of HiC data for all samples (paired-end)
├── publicData                      # Public HiC data that may be used for comparison
├── reference.files                 # reference files, mostly conda envs for tools used
│   ├── cooltools.env.yml
│   ├── distiller.env.yml
│   ├── TADLib.env.yml
│   └── README.md
├── results.NSC                     # HiC results produced by distiller pipeline and downstream analyses
│   ├── coolers_library
│   ├── fastqc
│   ├── mapped_parsed_sorted_chunks
│   └── pairs_library
├── sample.configs                  # Config files supplied to distiller to generate HiC matrices    
├── sample.metadata.tsv             # metadata for all HiC samples
└── scripts                         # all code
```
## Producing Results

### Notes

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

Use the `qc3C` tool in bam mode to profile quality metrics for our HiC samples.

```bash
$ ./scripts/matrix.utils.sh qc3c                             # run qc3C with specified params
        ./results.NSC/sample.QC/qc3C/                        # output dir
        results.NSC/mapped_parsed_sorted_chunks/16p.**/*.bam # bam files produced by distiller for each sample
# copy pastable
./scripts/matrix.utils.sh qc3c ./results.NSC/sample.QC/qc3C/ results.NSC/mapped_parsed_sorted_chunks/16p.**/*.bam
```
### Generate MultiQC reports

Distiller outputs multiple QC reports/files per sample than can each be aggregated into its own qc report with all samples together. This is how we generate all those reports

```bash
$ ./scripts/matrix.utils.sh multiqcs             # generate multiqc reports for different outputs from distiller,qc3C
        ./results.NSC/sample.QC/multiqc.reports/ # output dir
        ./results.NSC/                           # all multiqc data is under here in specific directories
# copy pastable
./scripts/matrix.utils.sh multiqcs ./results.NSC/sample.QC/multiqc.reports/ ./results.NSC/
```
