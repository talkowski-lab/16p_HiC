library(here)
library(glue)
###################################################
# Set Root Dirs
###################################################
if (grepl('/home/', BASE_DIR)) {
    SCRIPT_DIR <- file.path(BASE_DIR, '../remote.16p/scripts')
} else {
    SCRIPT_DIR <- here('scripts')
}
SAMPLE_METADATA_FILE <- file.path(BASE_DIR, 'HiC.16p.sample_metadata.tsv')
REF_DIR              <- file.path(BASE_DIR, 'reference.files')
RESULTS_DIR          <- file.path(BASE_DIR, 'results')

###################################################
# distiller-nf  
###################################################
# distiller-nf input
GENOME_REF_DIR          <- file.path(REF_DIR, 'genome.reference')
GENOME_REF_NAME         <- 'GRCh38_no_alt_analysis_set_GCA_000001405.15'
CHROMOSOME_SIZES_FILE   <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.chrom.sizes'))
BWA_INDEX_WILDCARD_PATH <- file.path(GENOME_REF_DIR, glue('{GENOME_REF_NAME}.fasta.*'))
# distiller-nf output 
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
###################################################
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
# HiCRep results
###################################################
HICREP_DIR          <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_DIR  <- file.path(HICREP_DIR, 'results')
HICREP_RESULTS_FILE <- file.path(HICREP_DIR, 'all.hicrep.scores.tsv')

###################################################
# Reproducing analysis from Weiner et al. 2022
###################################################
WEINER_REPLICATION_DIR       <- file.path(RESULTS_DIR, 'weiner.replication')
WEINER_REPLICATION_PLOTS_DIR <- file.path(WEINER_REPLICATION_DIR, 'plots')

###################################################
# TAD Annotations
###################################################
TAD_DIR                       <- file.path(RESULTS_DIR, 'TADs')
HITAD_TAD_RESULTS_FILE        <- file.path(TAD_DIR, 'all.hiTAD.TADs.tsv')
# HITAD_DI_RESULTS_FILE         <- file.path(TAD_DIR, 'all.hiTAD.DI.annotations.tsv')
HITAD_MOC_FILE                <- file.path(TAD_DIR, 'all.hiTAD.TAD.MoCs.tsv')
# COOLTOOLS_TAD_RESULTS_FILE    <- file.path(TAD_DIR, 'all.cooltools.TADs.tsv')
# COOLTOOLS_DI_RESULTS_FILE     <- file.path(TAD_DIR, 'all.cooltools.DI.annotations.tsv')
# COOLTOOLS_DI_RESULTS_FILE     <- file.path(TAD_DIR, 'all.cooltools.TAD.MoCs.tsv')
CONSENSUSTAD_TAD_RESULTS_FILE <- file.path(TAD_DIR, 'all.ConsensusTAD.TADs.tsv')
# CONSENSUSTAD_MOC_FILE         <- file.path(TAD_DIR, 'all.ConsensusTAD.DI.annotations.tsv')
CONSENSUSTAD_MOC_FILE         <- file.path(TAD_DIR, 'all.ConsensusTAD.TAD.MoCs.tsv')
# TADCompare results
TADCOMPARE_DIR                <- file.path(TAD_DIR, 'results_TADCompare')
TADCOMPARE_RESULTS_FILE       <- file.path(TAD_DIR, 'all.TADCompare.results.tsv')

###################################################
# Loop Annotations
###################################################
LOOPS_DIR                                  <- file.path(RESULTS_DIR, 'loops')
ALL_COOLTOOLS_LOOPS_RESULTS_FILE           <- file.path(LOOPS_DIR, 'all.cooltools.loops.tsv')
FILTERED_COOLTOOLS_LOOPS_RESULTS_FILE      <- file.path(LOOPS_DIR, 'filtered.cooltools.loops.tsv')
LOOPS_IDR2D_DIR                            <- file.path(LOOPS_DIR, 'results_IDR2D')
ALL_LOOPS_IDR2D_RESULTS_FILE               <- file.path(LOOPS_DIR, 'all.cooltools.IDR2D.results.tsv')
FILTERED_LOOPS_FILTERED_IDR2D_RESULTS_FILE <- file.path(LOOPS_DIR, 'filtered.cooltools.IDR2D.results.tsv')

###################################################
# Differential Contact results from multiHiCCompare
###################################################
MULTIHICCOMPARE_DIR                   <- file.path(RESULTS_DIR, 'multiHiCCompare')
ALL_MULTIHICCOMPARE_RESULTS_FILE      <- file.path(MULTIHICCOMPARE_DIR, 'all.multiHiCCompare.results.tsv')
FILTERED_MULTIHICCOMPARE_RESULTS_FILE <- file.path(MULTIHICCOMPARE_DIR, 'filtered.multiHiCCompare.results.tsv')
MULTIHICCOMPARE_SIG_RESULTS_FILE      <- file.path(MULTIHICCOMPARE_DIR, 'all.multiHiCCompare.n.results.tsv')

###################################################
# Compartment Annotations
###################################################
COMPARTMENTS_DIR         <- file.path(RESULTS_DIR, 'compartments')
COMPARTMENTS_RESULTS_DIR <- file.path(COMPARTMENTS_DIR, 'results')
# DCHIC_REF_DIR                  <- file.path(REF_DIR, 'dcHiC')
# COMPARTMENTS_PREPROCESSED_DIR  <- file.path(COMPARTMENTS_DIR, 'pre.processed.input')

###################################################
# gghic results
###################################################
GGHIC_DIR       <- file.path(RESULTS_DIR, 'gghic.plots')
GGHIC_PLOTS_DIR <- file.path(GGHIC_DIR, 'plots')

