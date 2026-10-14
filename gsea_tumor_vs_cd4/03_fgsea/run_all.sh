#!/bin/bash
# Run fgsea for every DE table in <de_dir>.
# Output: <out_dir>/<group>/set_<collection>.Rdata
#
# Usage: bash run_all.sh <de_dir> <out_dir>
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

for f in "$1"/*.csv; do
    g=$(basename "$f" .csv)
    Rscript "$HERE/run_fgsea.R" "$f" "$2/$g"
done
