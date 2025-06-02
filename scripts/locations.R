library(here)
# BASE_DIR <- here()
# Input Data
GENOME_REF_DIR                 <- '/data/talkowski/tools/ref/Hi_c_noalt'
CHROMOSOME_SIZES_FILE          <- file.path(GENOME_REF_DIR, 'hg38.reduced.chrom.sizes')
FASTQ_DIR                      <- file.path(BASE_DIR, 'fastq')
SAMPLE_METADATA_FILE           <- file.path(BASE_DIR, 'sample.metadata.tsv')
# distiller produced ou           puts
RESULTS_DIR                    <- file.path(BASE_DIR, 'results')
BAM_DIR                        <- file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
PAIRS_DIR                      <- file.path(RESULTS_DIR, 'pairs_library')
COOLERS_DIR                    <- file.path(RESULTS_DIR, 'coolers_library')
# Used for Sample QCing           
SAMPLE_QC_DIR                  <- file.path(RESULTS_DIR, 'sample.QC')
COVERAGE_DIR                   <- file.path(SAMPLE_QC_DIR, 'coverage')
# COVERAGE_DATA_FILE <-             file.path(COVERAGE_DIR, 'all.coverage.data.tsv')
# hicrep results
HICREP_DIR                     <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE            <- file.path(HICREP_DIR, 'all.hicrep.scores.tsv')
# Reproducing figure from Elise Robinson
ROBINSON_REPLICATION_DIR       <- file.path(RESULTS_DIR, 'Robinson.replication')
ROBINSON_REPLICATION_DATA_FILE <- file.path(ROBINSON_REPLICATION_DIR, 'replication.data.chr16.tsv')
# Differential Contact results from multiHiCCompare
SPARSE_MATRIX_DIR              <- file.path(RESULTS_DIR, 'sparse.matrices')
MULTIHICCOMPARE_DIR            <- file.path(RESULTS_DIR, 'multiHiCCompare')
# TAD Annotations
TAD_DIR                        <- file.path(RESULTS_DIR, 'TADs')
TAD_RESULTS_FILE               <- file.path(TAD_DIR, 'all.TAD.annotations.tsv')
