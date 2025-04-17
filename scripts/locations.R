library(glue)
library(here)
BASE_DIR <- here()
SAMPLE_METADATA_TSV <- glue('{BASE_DIR}/sample.metadata.tsv')
# Input Data
FASTQ_DIR <- 
    glue('{BASE_DIR}/fastq')
COOLERS_DIR <- 
    glue('{BASE_DIR}/results/coolers_library')
SPARSE_MATRIX_DIR <- 
    glue('{BASE_DIR}/results/sparse.matrices')
MULTIHICCOMPARE_DIR <- 
    glue('{BASE_DIR}/results/multiHiCCompare.results')
ELISE_RECREATION_DIR <-  
    glue('{BASE_DIR}/results/Elise.Recreation.results')
MATRIX_QC_DIR <-  
    glue('{BASE_DIR}/results/Matrix.QC')
HICREP_RESULTS_DIR <- 
    glue('{BASE_DIR}/results/hicrep.results')
TAD_ANNOTATIONS_DIR <- 
    glue("{BASE_DIR}/results/TADs")
# Factor levels for variaous metadata categories
EDITS <- 
    c(
      '16p',
      'WAPL',
      'NIBPL',
      'RAD21'
    )
GENOTYPES <- 
    c(
      'WT',
      'DEL',
      'DUP'
    )
CELLTYPES <- 
    c(
      'NSC',
      'iN'
    )
