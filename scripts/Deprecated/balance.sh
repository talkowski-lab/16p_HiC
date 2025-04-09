#!/bin/bash
set -euo pipefail

HIC_SAMPLES=${@}
for sample_file in ${HIC_SAMPLES[@]}; do
    sample_file="$(readlink -e ${sample_file})"
    for uri in $(cooler ls "$sample_file"); do
        cooler balance $uri
        echo $uri
    done
done
