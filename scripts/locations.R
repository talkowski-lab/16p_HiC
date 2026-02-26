library(here)
library(glue)
###################################################
# Script dir location
###################################################
if (grepl('/home/', BASE_DIR)) {
    SCRIPT_DIR <- file.path(BASE_DIR, '../remote.16p/scripts')
} else {
    SCRIPT_DIR <- here('scripts')
}

###################################################
# distiller-nf  
###################################################
SAMPLE_METADATA_FILE    <- file.path(BASE_DIR, 'HiC.16p.sample_metadata.tsv')
REF_DIR                 <- file.path(BASE_DIR, 'reference.files')
GENOME_REF_DIR          <- file.path(REF_DIR, 'genome.reference')
GENOME_REF_NAME         <- 'GRCh38_no_alt_analysis_set_GCA_000001405.15'
CHROMOSOME_SIZES_FILE   <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.chrom.sizes'))
BWA_INDEX_WILDCARD_PATH <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.fasta.*'))
# distiller-nf ouput 
RESULTS_DIR             <- file.path(BASE_DIR, 'results')
BAM_DIR                 <- file.path(RESULTS_DIR, 'mapped_parsed_sorted_chunks')
PAIRS_DIR               <- file.path(RESULTS_DIR, 'pairs_library')
COOLERS_DIR             <- file.path(RESULTS_DIR, 'coolers_library')

###################################################
# Functional genome annotations
###################################################
ENCODE_CCRE_DIR              <- '/data/talkowski/tools/ref/ENCODE/cCRE/v4'
ENCODE_CCRE_ANNOTATIONS_FILE <- file.path(REF_DIR, 'ENCODE.v4.cCRE.anontations.tsv')
CTCF_SITE_FILE               <- file.path(REF_DIR, 'CTCF.annotations.tsv')
GENE_DESERT_REGIONS_FILE     <- file.path(REF_DIR, 'hg38_gene_deserts.tsv')
GENOME_GTF_FILE              <- file.path(REF_DIR, 'gencode.v38.annotation.gtf.gz')
GENE_CONSTRAINTS_FILE        <- file.path(REF_DIR, 'gene.constraints.CNVR.tsv')
# FUNCTIONAL_ANNOTATIONS_DIR <- '/data/talkowski/xuefang/data/gnomad_V3/module08/step16_reannotate/noncoding_analyses/nc_elements/encode3'
# ABC_ANNOTATIONS_FILE       <- file.path(FUNCTIONAL_ANNOTATIONS_DIR, '../abc', 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz')

###################################################
# RNASeq
EXPRESSION_DATA_DIR     <- file.path(RESULTS_DIR, 'RNASeq', 'expression')
EXPRESSION_RESULTS_FILE <- file.path(RESULTS_DIR, 'RNASeq', 'all.expression.data.tsv')
DESEQ2_DATA_DIR         <- file.path(RESULTS_DIR, 'RNASeq', 'DESeq2')
DESEQ2_RESULTS_FILE     <- file.path(RESULTS_DIR, 'RNASeq', 'all.DESeq2.results.tsv')
TRADE_RESULTS_FILE      <- file.path(RESULTS_DIR, 'RNASeq', 'all.TRADE.results.tsv')

###################################################
# Sample QC Results
###################################################
SAMPLE_QC_DIR                   <- file.path(RESULTS_DIR, 'sample.QC')
COVERAGE_DIR                    <- file.path(SAMPLE_QC_DIR, 'coverage')
RESOLTION_COVERAGE_SUMAMRY_FILE <- file.path(SAMPLE_QC_DIR, 'resolution.coverage.summaries.tsv')
MIN_SAMPLE_RESOLUTION_FILE      <- file.path(SAMPLE_QC_DIR, 'minimum.viable.resolutions.tsv')

###################################################
# hicrep results
###################################################
HICREP_DIR                     <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_FILE            <- file.path(HICREP_DIR, 'all.hicrep.scores.tsv')

###################################################
# Reproducing figure from Elise Robinson
###################################################
ROBINSON_REPLICATION_DIR       <- file.path(RESULTS_DIR, 'robinson.replication')
ROBINSON_REPLICATION_DATA_FILE <- file.path(ROBINSON_REPLICATION_DIR, 'replication.data.chr16.tsv')

###################################################
# Differential Contact results from multiHiCCompare
###################################################
MULTIHICCOMPARE_DIR            <- file.path(RESULTS_DIR, 'multiHiCCompare')
MULTIHICCOMPARE_RESULTS_FILE   <- file.path(MULTIHICCOMPARE_DIR, 'multiHiCCompare.results.tsv')
MULTIHICCOMPARE_SIG_RESULTS_FILE <- file.path(MULTIHICCOMPARE_DIR, 'multiHiCCompare.n.results.tsv')

###################################################
# TAD Annotations
###################################################
TAD_DIR                         <- file.path(RESULTS_DIR, 'TADs')
HITAD_TAD_RESULTS_FILE          <- file.path(TAD_DIR, 'all.hiTAD.TADs.tsv')
# HITAD_DI_RESULTS_FILE          <- file.path(TAD_DIR, 'all.hiTAD.DI.annotations.tsv')
# HITAD_MOC_FILE                 <- file.path(TAD_DIR, 'all.hiTAD.TAD.MoCs.tsv')
COOLTOOLS_TAD_RESULTS_FILE      <- file.path(TAD_DIR, 'all.cooltools.TADs.tsv')
# COOLTOOLS_DI_RESULTS_FILE      <- file.path(TAD_DIR, 'all.cooltools.DI.annotations.tsv')
CONSENSUSTAD_TAD_RESULTS_FILE  <- file.path(TAD_DIR, 'all.ConsensusTAD.TADs.tsv')
# CONSENSUSTAD_MOC_FILE          <- file.path(TAD_DIR, 'all.ConsensusTAD.TAD.MoCs.tsv')
# TADCOMPARE_TAD_INPUT_FILE      <- file.path(TAD_DIR, 'all.TADCompare.TAD.inputs.tsv')
ALL_TAD_RESULTS_FILE            <- file.path(TAD_DIR, 'all.TADs.tsv')
TADCOMPARE_RESULTS_FILE         <- file.path(TAD_DIR, 'all.TADCompare.results.tsv')
ALL_TAD_SIMILARITY_RESULTS_FILE <- file.path(TAD_DIR, 'all.TAD.MoCs.tsv')

###################################################
# Compartment Annotations
###################################################
# DCHIC_REF_DIR                  <- file.path(REF_DIR, 'dcHiC')
# COMPARTMENTS_DIR               <- file.path(RESULTS_DIR, 'compartments')
# COMPARTMENTS_PREPROCESSED_DIR  <- file.path(COMPARTMENTS_DIR, 'pre.processed.input')
# COMPARTMENTS_RESULTS_DIR       <- file.path(COMPARTMENTS_DIR, 'results')

###################################################
# Loop Annotations
###################################################
LOOPS_DIR                     <- file.path(RESULTS_DIR, 'loops')
COOLTOOLS_LOOPS_RESULTS_FILE  <- file.path(LOOPS_DIR, 'all.cooltools.loops.tsv')
LOOPS_IDR2D_DIR               <- file.path(LOOPS_DIR, 'IDR2D')
LOOPS_IDR2D_RESULTS_FILE      <- file.path(LOOPS_DIR, 'all.cooltools.IDR2D.results.tsv')
LOOPS_IDR2D_GENE_MAPPING_FILE <- file.path(LOOPS_DIR, 'all.cooltools.IDR2D.and.genes.tsv')
# LOOP_SIMILARITY_RESULTS_FILE <- file.path(LOOPS_DIR, 'all.cooltools.loop.similarities.tsv')

###################################################
# gghic objects
###################################################
GGHIC_DIR <- file.path(RESULTS_DIR, 'gghic.plots')

