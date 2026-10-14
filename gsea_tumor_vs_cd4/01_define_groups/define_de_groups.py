"""Define the cell groups compared by differential expression (tumor cells vs normal CD4 T cells).

Reproduces the tumor / CD4 labels of the cell-annotation files used by the original MAST runs:
  CPDB/NOS_vs_Normal/normal_reanno.csv (= CPDB/AITL_vs_Normal/normal_reanno.csv) -> CD4 reference
  CPDB/AITL_vs_Normal/aitl_reanno.csv  (Aug-2025 version, with T_tumor_135)      -> TFHL (AITL) tumor
  CPDB/NOS_vs_Normal/nos_reanno.csv                                              -> PTCL-NOS tumor
  CPDB/EBV/ebv_reanno.csv                                                        -> EBV+ PTCL tumor
Only the metadata tables are needed (the expression matrix is not read here).

Usage:
  python define_de_groups.py --data-dir <dir containing Integration/ and Subtype/> [--out de_groups.csv]

Output: CSV with columns cell, group. group is 'CD4' or the tumor group name
(AITL_snXX, NOS_snXX, EBV_snXX); cells not used by DE are not listed.
"""
import argparse
import os

import pandas as pd

NORMAL_SAMPLES = ['sn04', 'sn08', 'sn09', 'sn10', 'sn15', 'sn16']
AITL_SAMPLES = ['sn01', 'sn02', 'sn03', 'sn05', 'sn06', 'sn07', 'sn12', 'sn13', 'sn14',
                'sn29', 'sn31', 'sn35', 'sn36', 'sn39', 'sn43', 'sn45', 'sn46', 'sn49']
NOS_SAMPLES = ['sn30', 'sn32', 'sn33', 'sn34', 'sn37', 'sn40', 'sn41', 'sn44', 'sn48']
EBV_SAMPLES = ['sn38', 'sn47']

# Major annotation of the integrated leiden_res0.7 clusters
MAJOR = {
    '0': 'T-cell', '1': 'B-cell', '2': 'T-cell', '3': 'plasma', '4': 'T-cell',
    '5': 'macrophage', '6': 'T-cell', '7': 'stromal', '8': 'NOS_sn34', '9': 'NOS_sn41',
    '10': 'NOS_sn37', '11': 'NOS_sn40', '12': 'NOS_sn48', '13': 'NOS_sn32',
    '14': 'endothelium', '15': 'NOS_sn44', '16': 'unknown16', '17': 'GCB',
    '18': 'monocyte', '19': 'unknown19', '20': 'pDC', '21': 'unknown21',
}
UNKNOWN_CLUSTERS = ['16', '19', '21']


def load_meta(data_dir):
    """Per-cell batch and major annotation (unknown clusters removed)."""
    meta = pd.read_csv(os.path.join(data_dir, 'Integration/meta.csv'), index_col=0)
    obs = pd.DataFrame(index=meta.index)
    obs['batch'] = [i[:4] for i in obs.index]
    obs['cluster'] = meta['leiden_res0.7'].astype(str).values
    obs = obs[~obs['cluster'].isin(UNKNOWN_CLUSTERS)].copy()
    obs['annotation'] = [MAJOR[c] for c in obs['cluster']]
    return obs


def load_subclusters(data_dir):
    rd = lambda p: pd.read_csv(os.path.join(data_dir, p), index_col=0)
    df_T = rd('Subtype/T-cell/meta.csv')
    df_l14 = rd('Subtype/leiden_14/meta.csv')
    df_l7 = rd('Subtype/leiden_7/meta.csv')
    df_mac = rd('Subtype/Macrophage/meta.csv')
    return {
        'T': df_T['leiden_res0.7'],
        'T_batch': df_T['batch'],
        'mast': df_l14.index[df_l14['leiden_res0.9'] == 10],
        'cDC': df_mac.index[df_mac['leiden_res0.1'] == 1],
        'to_remove': list(df_l14.index[df_l14['leiden_res0.9'].isin([4, 5, 7, 11])]) +
                     list(df_l7.index[df_l7['leiden_res0.8'].isin([2, 7, 10, 11])]),
    }


def relabel(obs, groups, sub):
    """Drop low-quality cells, then overwrite the annotation with sub-cluster labels.

    groups: list of (label, cell_index); later groups win on overlap.
    """
    obs = obs[~obs.index.isin(sub['to_remove'])].copy()
    dct = {}
    for label, cells in groups:
        dct.update(dict.fromkeys(cells, label))
    obs['reanno'] = [dct.get(c, a) for c, a in zip(obs.index, obs['annotation'])]
    return obs


def t_groups(sub, cd8_clusters, tumor_groups):
    t = sub['T']
    return ([('CD4', t.index[t == 0])] + tumor_groups +
            [('Treg', t.index[t == 3]),
             ('CD8', t.index[t.isin(cd8_clusters)]),
             ('mast', sub['mast']),
             ('cDC', sub['cDC'])])


def is_nos(labels):
    return labels.str.split('_').str[0] == 'NOS'


def normal_cd4(obs, sub):
    """CD4 reference: T sub-cluster 0 of the 6 reactive lymph nodes, pooled."""
    t = sub['T']
    normal = obs[obs['batch'].isin(NORMAL_SAMPLES)]
    normal = normal[~is_nos(normal['annotation'])]
    normal = relabel(normal, t_groups(sub, [1, 5, 12], [('TFH', t.index[t.isin([2, 4])])]), sub)
    return normal.index[normal['reanno'] == 'CD4']


def aitl_tumor(obs, sub):
    """TFHL tumor: T sub-clusters 2 and 4; in sn13 also T sub-cluster 5.

    The sn13 sub-cluster-5 cells (296) fall in an integrated cluster annotated NOS_sn37 and are
    kept, while all other NOS_* cells in AITL samples are dropped. This was added in Aug 2025
    and only changes sn13 (see README).
    """
    t, tb = sub['T'], sub['T_batch']
    sn13_tumor = t.index[(t == 5) & (tb == 'sn13')]
    aitl = obs[obs['batch'].isin(AITL_SAMPLES)]
    aitl = aitl[~is_nos(aitl['annotation']) | aitl.index.isin(sn13_tumor)]
    tumor = [('T_tumor', t.index[t.isin([2, 4])]), ('T_tumor', sn13_tumor)]
    aitl = relabel(aitl, t_groups(sub, [1, 12], tumor), sub)
    aitl = aitl[aitl['reanno'] == 'T_tumor']
    return 'AITL_' + aitl['batch']


def nos_tumor(obs, sub):
    """PTCL-NOS tumor: sample-specific integrated clusters (NOS_snXX), plus
    T sub-cluster 11 in sn33 and the Treg-like cells of sn30.

    As in the original analysis, a NOS_snXX group holds every NOS-sample cell in that cluster,
    including a few cells from other NOS samples (e.g. 9 of the 3,680 NOS_sn32 cells).
    """
    t, tb = sub['T'], sub['T_batch']
    nos = obs[obs['batch'].isin(NOS_SAMPLES)]
    nos = relabel(nos, t_groups(sub, [1, 12], [('TFH', t.index[t.isin([2, 4])])]), sub)
    sn33_tumor = t.index[(t == 11) & (tb == 'sn33')]
    nos.loc[nos.index.isin(sn33_tumor), 'reanno'] = 'NOS_sn33'
    nos.loc[(nos['batch'] == 'sn30') & (nos['reanno'] == 'Treg'), 'reanno'] = 'NOS_sn30'
    nos = nos[is_nos(nos['reanno'])]
    return nos['reanno']


def ebv_tumor(obs, sub):
    """EBV+ PTCL tumor: T sub-cluster 12 within sn38 and sn47."""
    t, tb = sub['T'], sub['T_batch']
    ebv = obs[obs['batch'].isin(EBV_SAMPLES)]
    ebv = relabel(ebv, t_groups(sub, [1, 12], [('T_tumor', t.index[t.isin([2, 4])])]), sub)
    tumor = t.index[(t == 12) & tb.isin(EBV_SAMPLES)]
    ebv = ebv[ebv.index.isin(tumor)]
    return 'EBV_' + ebv['batch']


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', required=True,
                    help='Folder containing Integration/ and Subtype/ (upstream snRNA-seq analysis)')
    ap.add_argument('--out', default='de_groups.csv')
    args = ap.parse_args()

    obs = load_meta(args.data_dir)
    sub = load_subclusters(args.data_dir)
    groups = pd.concat([
        pd.Series('CD4', index=normal_cd4(obs, sub)),
        aitl_tumor(obs, sub),
        nos_tumor(obs, sub),
        ebv_tumor(obs, sub),
    ])
    assert groups.index.is_unique
    groups.rename_axis('cell').rename('group').to_csv(args.out)
    print(groups.value_counts().sort_index().to_string())
