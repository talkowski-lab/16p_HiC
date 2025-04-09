# 16p_HiC
Analysis and code of HiC data for 16p editied cell lines

Current relevant file structure
```
./
├── distiller-nf                        # installed distiller nextflow pipeline
├── fastq                               # fastq files of HiC data for all samples (paired-end)
├── publicData                          # Public HiC data that may be used for comparison
├── reference.files                     # reference files, mostly conda envs for tools used
│   ├── cooltools.env.yml
│   ├── distiller.env.yml
│   ├── TADLib.env.yml
│   └── README.md
├── results.NSC                         # HiC results produced by distiller pipeline and downstream analyses
│   ├── coolers_library
│   ├── fastqc
│   ├── mapped_parsed_sorted_chunks
│   └── pairs_library
├── sample.configs                      # Config files supplied to distiller to generate HiC matrices    
├── sample.metadata.tsv                 # metadata for all HiC samples
└── scripts                             # all code
```
