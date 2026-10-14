#!/bin/bash
# Run MAST for the 8 comparisons: 4 TME cell types x {TFHL (AITL), PTCL-NOS}, each vs the same
# cell type in the pooled reactive lymph nodes.
# The original jobs ran on SLURM, one Rscript call per comparison.
#
# Usage: bash run_all.sh <raw_cellbender_DF_cat.h5ad> <de_groups.csv> <out_dir>
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

Rscript "$HERE/run_mast.R" "$1" "$2" "$3" \
    Macrophage-AITL Plasma-AITL B-cell-AITL cDC-AITL \
    Macrophage-NOS  Plasma-NOS  B-cell-NOS  cDC-NOS
