"""Generate a small simulated snRNA-seq dataset for the PTCL_atlas demo.

The real data are under controlled access in dbGaP (phs004552.v1.p1), so the
demo uses fully simulated counts instead of patient data. The simulation mimics
the structure the pipeline expects: several "samples" (batches) with a
sample-specific batch effect, and a handful of shared cell types.

Outputs (written to the current directory):
  raw_cellbender_DF_cat.h5ad  input for snRNA-seq_initial_clustering/run_scanorama.py
  raw_data.h5ad               input for snRNA-seq_subclustering/SubClustering.py
  random_select.pkl           down-sampling masks for SubClustering.py

Usage:
  python make_demo_data.py [--n-samples 12] [--cells-per-sample 300] [--n-subsamples 10]
"""
import argparse
import pickle

import anndata as ad
import numpy as np
import pandas as pd


def simulate(n_samples, cells_per_sample, n_genes, n_celltypes, seed):
    rng = np.random.default_rng(seed)
    # Gene names: a few mitochondrial genes so QC metrics behave as on real data
    genes = [f"MT-GENE{i}" for i in range(13)] + [f"GENE{i}" for i in range(n_genes - 13)]

    base = rng.gamma(shape=0.6, scale=1.0, size=n_genes)  # baseline expression
    # Each cell type up-regulates its own block of marker genes
    type_profiles = np.tile(base, (n_celltypes, 1))
    block = (n_genes - 13) // n_celltypes
    for t in range(n_celltypes):
        idx = 13 + t * block + rng.choice(block, size=60, replace=False)
        type_profiles[t, idx] *= rng.uniform(4, 10, size=60)

    mats, obs = [], []
    for s in range(n_samples):
        batch = f"sn{s + 1:02d}"
        batch_effect = rng.lognormal(mean=0.0, sigma=0.3, size=n_genes)
        # Unequal cell-type proportions per sample
        props = rng.dirichlet(np.ones(n_celltypes) * 2)
        types = rng.choice(n_celltypes, size=cells_per_sample, p=props)
        depth = rng.lognormal(mean=np.log(2000), sigma=0.4, size=cells_per_sample)
        mu = type_profiles[types] * batch_effect
        mu = mu / mu.sum(axis=1, keepdims=True) * depth[:, None]
        counts = rng.poisson(mu).astype(np.float32)
        mats.append(counts)
        # Cell barcodes follow the "<batch>_<barcode>" convention used by run_prepare.py
        obs.append(pd.DataFrame(
            {"true_celltype": [f"type{t}" for t in types]},
            index=[f"{batch}_CELL{c:05d}-1" for c in range(cells_per_sample)],
        ))

    X = np.vstack(mats)
    obs = pd.concat(obs)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    return adata


def make_random_select(obs_names, n_subsamples, fraction, seed):
    """Boolean masks in the format SubClustering.py reads: columns '<PERC>_<ORDER>'.

    '1.0_1' keeps every cell (the reference clustering); '0.9_1' ... '0.9_N'
    each keep a random 90% of cells.
    """
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(index=obs_names)
    df["1.0_1"] = True
    for order in range(1, n_subsamples + 1):
        df[f"{fraction}_{order}"] = rng.random(len(obs_names)) < fraction
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-samples", type=int, default=12,
                   help="number of simulated samples; must be >= 10 because the pipeline keeps "
                        "genes that are highly variable in at least 10 samples (default 12)")
    p.add_argument("--cells-per-sample", type=int, default=300)
    p.add_argument("--n-genes", type=int, default=3000)
    p.add_argument("--n-celltypes", type=int, default=5)
    p.add_argument("--n-subsamples", type=int, default=10,
                   help="number of 90%% down-samplings for the stability analysis (paper used 100)")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    if args.n_samples < 10:
        p.error("--n-samples must be at least 10")

    adata = simulate(args.n_samples, args.cells_per_sample, args.n_genes, args.n_celltypes, args.seed)
    adata.write("raw_cellbender_DF_cat.h5ad")

    # SubClustering.py expects raw counts with a 'batch' column
    sub = adata.copy()
    sub.obs["batch"] = [i[:4] for i in sub.obs_names]
    sub.write("raw_data.h5ad")

    masks = make_random_select(sub.obs_names, args.n_subsamples, 0.9, args.seed)
    with open("random_select.pkl", "wb") as fh:
        pickle.dump(masks, fh)

    print(f"Wrote {adata.n_obs} cells x {adata.n_vars} genes from {args.n_samples} simulated samples")
    print("Files: raw_cellbender_DF_cat.h5ad, raw_data.h5ad, random_select.pkl")


if __name__ == "__main__":
    main()
