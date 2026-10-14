#!/bin/bash
# Make the ligand-activity scatter plot (and table) of each receiver cell type.
#
# Usage: bash run_all.sh <activity_dir> <out_dir>
#   activity_dir: output of 02_ligand_activity/run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$2"

# receiver : file-name tag of the original PDFs
for x in macrophage:macro endothelium:endothelium stromal:stromal; do
    receiver=${x%%:*}
    tag=${x##*:}
    Rscript "$HERE/plot_ligand_activity_nq.R" "$1" "$receiver" \
        "$2/plot_ligand_activity_nq_${tag}_fixed_goi.pdf" "$2/aupr_corrected_${receiver}.csv"
done
