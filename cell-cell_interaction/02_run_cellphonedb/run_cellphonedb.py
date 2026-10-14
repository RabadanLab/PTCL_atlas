"""Run CellPhoneDB statistical analysis on one dataset.

Replaces the per-folder run_cellphonedb.py scripts and ChangeName.sh of the original analysis.
Output files are renamed to <out_dir>/{means,pvalues,significant_means,deconvoluted,...}.txt

Usage:
  python run_cellphonedb.py --cpdb cellphonedb.zip --counts aitl_sn01.h5ad \
      --meta aitl_reanno_sn01.csv --out-dir cpdb_output/AITL/sn01
"""
import argparse
import glob
import os

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

ap = argparse.ArgumentParser()
ap.add_argument('--cpdb', required=True, help='CellPhoneDB database (cellphonedb-data v4.1.0 cellphonedb.zip)')
ap.add_argument('--counts', required=True, help='log-normalised counts (.h5ad), gene symbols as var_names')
ap.add_argument('--meta', required=True, help='CSV: cell barcode, cell type')
ap.add_argument('--out-dir', required=True)
ap.add_argument('--threads', type=int, default=4)
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

# Package defaults were used in the original analysis (stated explicitly here);
# debug_seed=-1 means the permutations are not seeded.
cpdb_statistical_analysis_method.call(
    cpdb_file_path=args.cpdb,
    meta_file_path=args.meta,
    counts_file_path=args.counts,
    counts_data='hgnc_symbol',
    output_path=args.out_dir,
    iterations=1000,
    threshold=0.1,
    threads=args.threads,
    result_precision=3,
    debug_seed=-1,
)

# statistical_analysis_<name>_<timestamp>.txt -> <name>.txt
for name in ['deconvoluted', 'deconvoluted_percents', 'means', 'pvalues', 'significant_means']:
    for f in glob.glob(os.path.join(args.out_dir, f'statistical_analysis_{name}_[0-9]*.txt')):
        os.replace(f, os.path.join(args.out_dir, f'{name}.txt'))
