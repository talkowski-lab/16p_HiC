#!/bin/bash

SIZES='/data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes'
REFERENCE='/data/talkowski/tools/ref/Hi-c_ref/hg38.fasta'
DIGESTED_REFERENCE="hg38.${ENZYME}.restricted.bed"
cooler digest $SIZES $REFERENCE $ENZYME > $DIGESTED_REFERENCE
pairtools restrict -f $DIGESTED_REFERENCE $pairs_file -o ${output_file}
