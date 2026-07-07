import scanpy as sc
import pandas as pd
import numpy as np
import scanorama
import math
from matplotlib import pyplot as plt
import sys
import random
from datetime import datetime
import pickle

#PERC=float(sys.argv[1])
PERC=sys.argv[1]
ORDER=sys.argv[2]


# Read in data
adata = sc.read('raw_data.h5ad')
df=pickle.load(open("random_select.pkl","rb"))


# Subselect
print("Random selection")
#RANDOM=np.random.binomial(1, PERC, adata.shape[0])>0
adata=adata[df.loc[:,f'{PERC}_{ORDER}'],:]
print(adata.shape)


# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Detect and select highly variable genes
sc.pp.highly_variable_genes(adata, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())

for I in range(len(set(adata.obs['batch']))+1):
    print(str(I)+"\t"+str(np.sum(adata.var.highly_variable_nbatches>=I)))

var_select = adata.var.highly_variable_nbatches >= 10
var_genes = var_select.index[var_select]
var_genes=[item for item in var_genes if not item.startswith('MT-')]
print(len(var_genes))


# Split data by batch
batches = np.unique(adata.obs['batch'])
alldata = {}
for BATCH in batches:
    alldata[BATCH] = adata[adata.obs['batch'] == BATCH,]


# Subset highly variable genes
alldata2 = dict()
for BATCH in alldata.keys():
    print(BATCH)
    alldata2[BATCH] = alldata[BATCH][:,var_genes]


#convert to list of AnnData objects
adatas = list(alldata2.values())


# Run Scanorama
scanorama.integrate_scanpy(adatas,dimred=50)

scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
all_s = np.concatenate(scanorama_int)
print(all_s.shape)


# add to the AnnData object, create a new object first
adata.obsm["Scanorama"] = all_s



# Calculate neighbors
sc.pp.neighbors(adata, n_pcs=30, use_rep="Scanorama", n_neighbors=50)
sc.tl.umap(adata)


# Leiden clustering
res_array = ['0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5']
df=pd.DataFrame(index=adata.obs.index)

for res in res_array:
    print(res)
    sc.tl.leiden(adata, resolution=float(res), key_added=f'leiden_res{res}', random_state = 42)
    df[f'res{res}']=adata.obs[f'leiden_res{res}']


# Output the results
df.to_csv(f'./output/subtypes_{PERC}_{ORDER}.csv')

if float(PERC)==1.0:
    adata.write(f'./output/integrated_{PERC}_{ORDER}.h5ad')
