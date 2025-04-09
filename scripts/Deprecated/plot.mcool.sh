#!/bin/bash
OUTPUT_DIR="/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p/results/Matrix.QC"
region="$1"
resolution="$2"
hic_samples=${@:3}
for sample_file in ${hic_samples[@]}; do
    sample_file="$(readlink -e ${sample_file})"
    sample_ID="$(basename "$sample_file")"
    sample_ID="${sample_ID%%.*.mcool}"
    # sample_ID="${sample_ID%%.hic}"
    output_file="${OUTPUT_DIR}/${sample_ID}.${resolution}.${region}.heatmap.pdf"
        # --weight-field cis.ice -vmin -1 -vmax 1 \
        # --weight-field weight -vmin 0 -vmax 1000 \
        # --weight-field cooler.balance -l \
    fancplot $region \
        -vv \
        --output ${output_file} \
        -p triangular \
        -u -l \
        --title "${sample_ID} Contacts on ${region%:*}" \
        ${sample_file}@${resolution}
    # fancplot chr16:0mb-96mb --output ./WT_VS_DEL.Merged.100K.heatmap.pdf -p split --title 'Merged WT vs Merged DEL Contacts chr16' ../coolers_library/16pWTNSCHIC_S5S6.hg38.mapq_30.1000.mcool::resolutions/100000 ../coolers_library/16pDELNSCHIC_S1S2.hg38.mapq_30.1000.mcool::resolutions/100000
done
