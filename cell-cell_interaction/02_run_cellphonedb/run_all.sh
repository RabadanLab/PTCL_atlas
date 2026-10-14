#!/bin/bash
# Run CellPhoneDB for the normal reference (pooled), each AITL/TFHL sample and each PTCL-NOS sample.
# Original jobs were submitted on SLURM with: --time=48:00:00 --mem=100G --cpus-per-task=10
#
# Usage: bash run_all.sh <cpdb_input_dir> <cpdb_output_dir> <path/to/cellphonedb.zip>
set -euo pipefail

IN=${1:-cpdb_input}
OUT=${2:-cpdb_output}
DB=${3:-cellphonedb.zip}
RUN="python3 $(dirname "$0")/run_cellphonedb.py --cpdb $DB"

# Normal reference (6 reactive lymph nodes pooled): provides the CD4 and TFH columns
$RUN --counts $IN/normal/normal.h5ad --meta $IN/normal/normal_reanno.csv --out-dir $OUT/normal

for S in sn01 sn02 sn03 sn05 sn06 sn07 sn12 sn13 sn14 sn29 sn31 sn35 sn36 sn39 sn43 sn45 sn46 sn49; do
    $RUN --counts $IN/AITL/aitl_$S.h5ad --meta $IN/AITL/aitl_reanno_$S.csv --out-dir $OUT/AITL/$S
done

for S in sn30 sn32 sn33 sn34 sn37 sn40 sn41 sn44 sn48; do
    $RUN --counts $IN/NOS/nos_$S.h5ad --meta $IN/NOS/nos_reanno_$S.csv --out-dir $OUT/NOS/$S
done
