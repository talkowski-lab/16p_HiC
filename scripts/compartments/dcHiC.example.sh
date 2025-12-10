COOLER_DIR="$(pwd)/results/coolers_library"

cool_file1="${COOLER_DIR}/16p.NSC.WT.Merged.Merged/16p.NSC.WT.Merged.Merged.hg38.mapq_30.1000.mcool"
prefix1="16p.NSC.WT"
cool_file2="${COOLER_DIR}/16p.NSC.DEL.Merged.Merged/16p.NSC.DEL.Merged.Merged.hg38.mapq_30.1000.mcool"
prefix2="16p.NSC.DEL"

./scripts/compartments/run.compartments.dcHiC.sh -r 100000 -t 1 -P "${prefix1}" -C "${cool_file1}" -p "${prefix2}" -c "${cool_file2}"
