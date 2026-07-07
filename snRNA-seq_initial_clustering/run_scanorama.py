import scanpy as sc
import pandas as pd
import numpy as np
import scanorama


# Read in data
adata = sc.read('raw_cellbender_DF_cat.h5ad')
adata.obs['batch']=[item[:4] for item in adata.obs.index]
adata.obs['batch'].isin(['sn38','sn47'])


# Calculate QC metrics
adata.var['mt'] = np.asarray([item.startswith('MT-') for item in adata.var.index])
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Detect and select highly variable genes
sc.pp.highly_variable_genes(adata, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())

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
adata_sc = adata.copy()
adata_sc.obsm["Scanorama"] = all_s


adata_sc.write('integrated.h5ad')
