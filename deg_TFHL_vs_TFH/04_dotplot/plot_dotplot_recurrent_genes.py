"""Dot plot of the recurrent TFHL-vs-TFH genes: normal TFH cells and TFHL tumor cells, one column per sample.
Concise version of ntbk_viz_TFHL_vs_TFH_genes.ipynb of the original analysis (the call writing
figures/dotplot_recurrent_genes_TFHL_vs_TFH.pdf).

The original notebook read the cell labels and log-normalised expression from
AnnotationLin/integrated_processed_annotation_Lin_Tcell.h5ad (built by AnnotationLin/ntbk_annotation_wendy.ipynb).
This script rebuilds the same cells and values from the upstream inputs:
  * cells: T sub-clusters 2 and 4 (in sn13 also T sub-cluster 5) of the 6 reactive lymph nodes (TFH)
    and of the 18 TFHL samples (tumor), without the low-quality sub-clusters and integrated clusters 16, 19.
    Unlike the MAST groups (01_define_groups), cells in the NOS_* integrated clusters are kept
    (151 more cells; see README).
  * expression: normalize_total (1e4 counts per cell) + log1p of the raw counts.

Usage:
  python plot_dotplot_recurrent_genes.py --data-dir <dir containing Integration/ and Subtype/>
                                         [--out dotplot_recurrent_genes_TFHL_vs_TFH.pdf]
"""
import argparse
import os

import anndata as ad
import pandas as pd
import scanpy as sc

NORMAL_SAMPLES = ['sn04', 'sn08', 'sn09', 'sn10', 'sn15', 'sn16']
# TFHL samples in the column order of the figure
AITL_SAMPLES = ['sn03', 'sn07', 'sn01', 'sn12', 'sn14', 'sn43', 'sn29',
                'sn06', 'sn39', 'sn02', 'sn05', 'sn13', 'sn49', 'sn45',
                'sn31', 'sn46', 'sn35', 'sn36']
BATCH_ORDER = NORMAL_SAMPLES + AITL_SAMPLES


def tfh_and_tumor_cells(data_dir):
    rd = lambda p: pd.read_csv(os.path.join(data_dir, p), index_col=0)
    cluster = rd('Integration/meta.csv')['leiden_res0.7'].astype(str)
    df_T = rd('Subtype/T-cell/meta.csv')
    df_l14 = rd('Subtype/leiden_14/meta.csv')
    df_l7 = rd('Subtype/leiden_7/meta.csv')
    to_remove = (set(df_l14.index[df_l14['leiden_res0.9'].isin([4, 5, 7, 11])]) |
                 set(df_l7.index[df_l7['leiden_res0.8'].isin([2, 7, 10, 11])]))

    t, batch = df_T['leiden_res0.7'], df_T.index.str[:4]
    keep = ((t.isin([2, 4]) | ((t == 5) & (batch == 'sn13'))) &
            batch.isin(BATCH_ORDER) &
            ~df_T.index.isin(to_remove) &
            ~cluster.reindex(df_T.index).isin(['16', '19']))
    return df_T.index[keep]


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', required=True,
                    help='Folder containing Integration/ and Subtype/ (upstream snRNA-seq analysis)')
    ap.add_argument('--out', default='dotplot_recurrent_genes_TFHL_vs_TFH.pdf')
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    cells = tfh_and_tumor_cells(args.data_dir)
    adata = ad.read_h5ad(os.path.join(args.data_dir, 'Integration/raw_cellbender_DF_cat.h5ad'), backed='r')
    adata = adata[adata.obs.index.isin(cells), :].to_memory()
    adata.obs['batch'] = pd.Categorical([i[:4] for i in adata.obs.index], categories=BATCH_ORDER)
    print(adata.obs['batch'].value_counts().reindex(BATCH_ORDER).to_string())

    # normalisation is per cell, so normalising after subsetting gives the original values
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # genes grouped by the number of samples in which they are significant
    selected = pd.read_csv(os.path.join(here, 'selected_genes.txt'), sep='\t', comment='#')
    var_names = {str(n): list(g['gene']) for n, g in selected.groupby('num_samples', sort=False)}

    dp = sc.pl.dotplot(adata, var_names=var_names, groupby='batch', swap_axes=True,
                       categories_order=BATCH_ORDER, standard_scale='var', cmap='OrRd',
                       return_fig=True)
    dp.savefig(args.out, bbox_inches='tight')
