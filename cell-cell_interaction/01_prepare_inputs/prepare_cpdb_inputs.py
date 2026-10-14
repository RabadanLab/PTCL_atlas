"""Prepare CellPhoneDB inputs (log-normalised counts + cell-type labels).

Consolidates the original notebooks:
  NOS_vs_Normal/ntbk_prepare_normal.ipynb   -> normal (reactive LN) reference, one run
  AITL_vs_Normal/ntbk_prepare_aitl.ipynb    -> AITL / TFHL (Mar-2024 version, see README)
  AITL_per_sample/ntbk_prepare_aitl.ipynb   -> split AITL per sample
  NOS_vs_Normal/ntbk_prepare_nos.ipynb      -> PTCL-NOS
  NOS_per_sample/ntbk_prepare_nos.ipynb     -> split NOS per sample

Usage:
  python prepare_cpdb_inputs.py --data-dir <dir containing Integration/ and Subtype/> --out-dir <dir>

Outputs (in --out-dir):
  normal/normal.h5ad, normal/normal_reanno.csv
  AITL/aitl_<sample>.h5ad, AITL/aitl_reanno_<sample>.csv
  NOS/nos_<sample>.h5ad,   NOS/nos_reanno_<sample>.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc

NORMAL_SAMPLES = ['sn04', 'sn08', 'sn09', 'sn10', 'sn15', 'sn16']
AITL_SAMPLES = ['sn01', 'sn02', 'sn03', 'sn05', 'sn06', 'sn07', 'sn12', 'sn13', 'sn14',
                'sn29', 'sn31', 'sn35', 'sn36', 'sn39', 'sn43', 'sn45', 'sn46', 'sn49']
NOS_SAMPLES = ['sn30', 'sn32', 'sn33', 'sn34', 'sn37', 'sn40', 'sn41', 'sn44', 'sn48']

# Major annotation of the integrated leiden_res0.7 clusters
MAJOR = {
    '0': 'T-cell', '1': 'B-cell', '2': 'T-cell', '3': 'plasma', '4': 'T-cell',
    '5': 'macrophage', '6': 'T-cell', '7': 'stromal', '8': 'NOS_sn34', '9': 'NOS_sn41',
    '10': 'NOS_sn37', '11': 'NOS_sn40', '12': 'NOS_sn48', '13': 'NOS_sn32',
    '14': 'endothelium', '15': 'NOS_sn44', '16': 'unknown16', '17': 'GCB',
    '18': 'monocyte', '19': 'unknown19', '20': 'pDC', '21': 'unknown21',
}
UNKNOWN_CLUSTERS = ['16', '19', '21']


def load_data(data_dir):
    """Raw counts -> normalize_total(1e4) -> log1p, with major annotation."""
    adata = sc.read(os.path.join(data_dir, 'Integration/raw_cellbender_DF_cat.h5ad'))
    adata.obs['batch'] = [i[:4] for i in adata.obs.index]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    meta = pd.read_csv(os.path.join(data_dir, 'Integration/meta.csv'), index_col=0)
    assert (adata.obs.index == meta.index).all()
    adata.obs['cluster'] = meta['leiden_res0.7'].astype(str).values
    adata = adata[~adata.obs['cluster'].isin(UNKNOWN_CLUSTERS), :].copy()
    adata.obs['annotation'] = [MAJOR[c] for c in adata.obs['cluster']]
    return adata


def load_subclusters(data_dir):
    """Sub-clustering results used to refine labels and remove low-quality cells."""
    rd = lambda p: pd.read_csv(os.path.join(data_dir, p), index_col=0)
    df_T = rd('Subtype/T-cell/meta.csv')
    df_l14 = rd('Subtype/leiden_14/meta.csv')
    df_l7 = rd('Subtype/leiden_7/meta.csv')
    df_mac = rd('Subtype/Macrophage/meta.csv')
    t = df_T['leiden_res0.7']
    return {
        'T': t,
        'T_batch': df_T['batch'],
        'mast': df_l14.index[df_l14['leiden_res0.9'] == 10],
        'cDC': df_mac.index[df_mac['leiden_res0.1'] == 1],
        'to_remove': list(df_l14.index[df_l14['leiden_res0.9'].isin([4, 5, 7, 11])]) +
                     list(df_l7.index[df_l7['leiden_res0.8'].isin([2, 7, 10, 11])]),
    }


def relabel(adata, groups, sub):
    """Drop low-quality cells, then overwrite 'annotation' with sub-cluster labels.

    groups: list of (label, cell_index); later groups win on overlap.
    """
    adata = adata[~adata.obs.index.isin(sub['to_remove']), :].copy()
    dct = {}
    for label, cells in groups:
        dct.update(dict.fromkeys(cells, label))
    adata.obs['reanno'] = [dct.get(c, a) for c, a in zip(adata.obs.index, adata.obs['annotation'])]
    return adata


def t_groups(sub, cd8_clusters, tumor_labels):
    t = sub['T']
    groups = [('CD4', t.index[t == 0])]
    groups += [(lab, t.index[t.isin(cl)]) for lab, cl in tumor_labels]
    groups += [('Treg', t.index[t == 3]),
               ('CD8', t.index[t.isin(cd8_clusters)]),
               ('mast', sub['mast']),
               ('cDC', sub['cDC'])]
    return groups


def write(adata, out_dir, name, tag=''):
    """Write <name><tag>.h5ad and <name>_reanno<tag>.csv (CellPhoneDB counts + meta)."""
    os.makedirs(out_dir, exist_ok=True)
    adata.write(os.path.join(out_dir, f'{name}{tag}.h5ad'))
    pd.DataFrame(adata.obs['reanno']).to_csv(os.path.join(out_dir, f'{name}_reanno{tag}.csv'))


def write_per_sample(adata, out_dir, name):
    for s in sorted(set(adata.obs['batch'])):
        write(adata[adata.obs['batch'] == s, :], out_dir, name, tag=f'_{s}')


def is_nos(labels):
    return np.asarray([l.split('_')[0] == 'NOS' for l in labels])


def prepare_normal(adata, sub, out_dir):
    normal = adata[adata.obs['batch'].isin(NORMAL_SAMPLES), :]
    normal = normal[~is_nos(normal.obs['annotation']), :]
    groups = t_groups(sub, cd8_clusters=[1, 5, 12], tumor_labels=[('TFH', [2, 4])])
    normal = relabel(normal, groups, sub)
    # Remaining unassigned 'T-cell' (T sub-clusters 10/11) are removed
    normal = normal[normal.obs['reanno'] != 'T-cell', :]
    write(normal, os.path.join(out_dir, 'normal'), 'normal')


def prepare_aitl(adata, sub, out_dir):
    aitl = adata[adata.obs['batch'].isin(AITL_SAMPLES), :]
    aitl = aitl[~is_nos(aitl.obs['annotation']), :]
    groups = t_groups(sub, cd8_clusters=[1, 5, 12],
                      tumor_labels=[('T_tumor', [2, 4])])
    aitl = relabel(aitl, groups, sub)
    aitl = aitl[aitl.obs['reanno'] != 'T-cell', :]
    write_per_sample(aitl, os.path.join(out_dir, 'AITL'), 'aitl')


def prepare_nos(adata, sub, out_dir):
    nos = adata[adata.obs['batch'].isin(NOS_SAMPLES), :]
    groups = t_groups(sub, cd8_clusters=[1, 12], tumor_labels=[('TFH', [2, 4])])
    nos = relabel(nos, groups, sub)
    t, tb = sub['T'], sub['T_batch']
    # sn33 tumor cells sit in T sub-cluster 11
    sn33_tumor = t.index[(t == 11) & (tb == 'sn33')]
    nos.obs.loc[nos.obs.index.isin(sn33_tumor), 'reanno'] = 'NOS_sn33'
    nos = nos[nos.obs['reanno'] != 'T-cell', :].copy()
    # In sn30 the tumor cells cluster with Treg
    nos.obs.loc[(nos.obs['batch'] == 'sn30') & (nos.obs['reanno'] == 'Treg'), 'reanno'] = 'NOS_sn30'
    write_per_sample(nos, os.path.join(out_dir, 'NOS'), 'nos')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', required=True,
                    help='Folder containing Integration/ and Subtype/ (upstream snRNA-seq analysis)')
    ap.add_argument('--out-dir', default='cpdb_input')
    args = ap.parse_args()

    adata = load_data(args.data_dir)
    sub = load_subclusters(args.data_dir)
    prepare_normal(adata, sub, args.out_dir)
    prepare_aitl(adata, sub, args.out_dir)
    prepare_nos(adata, sub, args.out_dir)
