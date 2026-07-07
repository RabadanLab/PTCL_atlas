import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import scanorama
import matplotlib.pyplot as plt
import anndata

import scipy.cluster
import scipy.cluster.hierarchy
import math
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
import scipy.stats as ss
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import fdrcorrection
from gap_statistic import OptimalK
from sklearn.cluster import KMeans
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist


adata=sc.read('integrated.h5ad')


sc.pp.neighbors(adata, n_pcs=30, use_rep="Scanorama", n_neighbors=50)
sc.tl.umap(adata)


Y=np.asarray(adata.obsm['Scanorama'][:,range(30)])
res_array = ['0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5']

for res in res_array:
    print(res)
    sc.tl.leiden(adata, resolution=float(res), key_added=f'leiden_res{res}', random_state = 42)

sil = [silhouette_score(Y, adata.obs[f'leiden_res{r}']) for r in res_array]
ks = [adata.obs[f'leiden_res{r}'].astype(int).max() + 1 for r in res_array]

print(sil)

plt.plot(res_array, sil, label='leiden')
for r, s, k in zip(res_array, sil, ks):
    plt.text(r, s*1.01, k, fontdict={'fontsize': 8})
plt.legend()
plt.title(f'n_neighbors: 50')
plt.show()


adata.write('integrated_processed.h5ad')
