#!/bin/bash
# Run MAST for all 29 tumor samples (18 TFHL, 9 PTCL-NOS, 2 EBV+ PTCL) vs pooled normal CD4 T cells.
# The original jobs ran on SLURM (24 h, 50 GB RAM, 5 CPUs), one Rscript call per sample.
#
# Usage: bash run_all.sh <raw_cellbender_DF_cat.h5ad> <de_groups.csv> <out_dir>
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

Rscript "$HERE/run_mast.R" "$1" "$2" "$3" \
    AITL_sn01 AITL_sn02 AITL_sn03 AITL_sn05 AITL_sn06 AITL_sn07 AITL_sn12 AITL_sn13 AITL_sn14 \
    AITL_sn29 AITL_sn31 AITL_sn35 AITL_sn36 AITL_sn39 AITL_sn43 AITL_sn45 AITL_sn46 AITL_sn49 \
    NOS_sn30 NOS_sn32 NOS_sn33 NOS_sn34 NOS_sn37 NOS_sn40 NOS_sn41 NOS_sn44 NOS_sn48 \
    EBV_sn38 EBV_sn47
