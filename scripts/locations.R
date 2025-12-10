library(here)
library(glue)
###################################################
# Script dir location
###################################################
if (grepl('/home/', BASE_DIR)) {
    SCRIPT_DIR <- here('../remote.16p/scripts')
} else {
    SCRIPT_DIR <- here('scripts')
}

###################################################
# distiller-nf 
###################################################
SAMPLE_METADATA_FILE    <- file.path(BASE_DIR, 'HiC.16p.sample_metadata.tsv')
# GENOME_REF_DIR          <- '/data/talkowski/tools/ref/Hi_c_noalt'
GENOME_REF_DIR          <- file.path(BASE_DIR, 'reference.files', 'genome.reference')
GENOME_REF_NAME         <- 'GRCh38_no_alt_analysis_set_GCA_000001405.15'
CHROMOSOME_SIZES_FILE   <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.chrom.sizes'))
BWA_INDEX_WILDCARD_PATH <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.fasta.*'))
# distiller-nf ouput 
RESULTS_DIR             <- file.path(BASE_DIR, 'results')
BAM_DIR                 <- file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
PAIRS_DIR               <- file.path(RESULTS_DIR, 'pairs_library')
COOLERS_DIR             <- file.path(RESULTS_DIR, 'coolers_library')

###################################################
# Sample QC Results
###################################################
SAMPLE_QC_DIR                   <- file.path(RESULTS_DIR, 'sample.QC')
COVERAGE_DIR                    <- file.path(SAMPLE_QC_DIR, 'coverage')
RESOLTION_COVERAGE_SUMAMRY_FILE <- file.path(SAMPLE_QC_DIR, 'resolution.coverage.summaries.tsv')
MIN_SAMPLE_RESOLUTION_FILE      <- file.path(SAMPLE_QC_DIR, 'minimum.viable.resolutions.data.tsv')

###################################################
# hicrep results
###################################################
HICREP_DIR                     <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE            <- file.path(HICREP_DIR, 'all.hicrep.scores.tsv')

###################################################
# Reproducing figure from Elise Robinson
###################################################
ROBINSON_REPLICATION_DIR       <- file.path(RESULTS_DIR, 'Robinson.replication')
ROBINSON_REPLICATION_DATA_FILE <- file.path(ROBINSON_REPLICATION_DIR, 'replication.data.chr16.tsv')

###################################################
# Differential Contact results from multiHiCCompare
###################################################
MULTIHICCOMPARE_DIR            <- file.path(RESULTS_DIR, 'multiHiCCompare')
MULTIHICCOMPARE_RESULTS_FILE   <- file.path(MULTIHICCOMPARE_DIR, 'multiHiCCompare.results.tsv')

###################################################
# TAD Annotations
###################################################
TAD_DIR                        <- file.path(RESULTS_DIR, 'TADs')
HITAD_TAD_RESULTS_FILE         <- file.path(TAD_DIR, 'all.hiTAD.TAD.annotations.tsv')
HITAD_DI_RESULTS_FILE          <- file.path(TAD_DIR, 'all.hiTAD.DI.annotations.tsv')
HITAD_MOC_FILE                 <- file.path(TAD_DIR, 'all.hiTAD.TAD.MoCs.tsv')

###################################################
# Functional genome annotations
###################################################
FUNCTIONAL_ANNOTATIONS_DIR     <- '/data/talkowski/xuefang/data/gnomad_V3/module08/step16_reannotate/noncoding_analyses/nc_elements'
FUNCTIONAL_ANNOTATION_FILES    <- list.files(file.path(FUNCTIONAL_ANNOTATIONS_DIR, 'encode3'))
ABC_ANNOTATIONS_FILE           <- file.path(FUNCTIONAL_ANNOTATIONS_DIR, 'abc', 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz')

