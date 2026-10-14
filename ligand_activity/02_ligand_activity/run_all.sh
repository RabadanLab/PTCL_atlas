#!/bin/bash
# Run NicheNet ligand activity for the 6 runs: {M1, M2} ecosystem gene sets x
# {macrophage, endothelium, stromal} receivers.
# The original analysis ran each in its own notebook.
#
# Usage: bash run_all.sh <raw_cellbender_DF_cat.h5ad> <meta_annotation_Lin.csv> <nichenet_dir> <out_dir>
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

Rscript "$HERE/run_nichenet.R" "$1" "$2" "$HERE/../01_genesets/genesets_of_interest.csv" "$3" "$4" \
    M1_to_macrophage  M2_to_macrophage \
    M1_to_endothelium M2_to_endothelium \
    M1_to_stromal     M2_to_stromal
