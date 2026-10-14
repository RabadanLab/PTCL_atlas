# Tumor cells vs normal CD4 T cells: differential expression and GSEA (snRNA-seq)

This folder contains the code that produces the heatmap **`plot_Fig3A_v3.pdf`** (Fig. 3A).
For each tumor sample, the tumor cells are compared with normal CD4 T cells from reactive lymph
nodes. The comparison is a differential expression test (MAST) followed by pre-ranked GSEA (fgsea).
The heatmap shows, for 29 selected gene sets, the signed significance
`-log10(padj) × sign(NES)` in each sample.

| Column group | Tumor cells | Samples |
|---|---|---|
| `AITL_*` | tumor cells of nodal TFH lymphoma (AITL / TFHL) | 18 |
| `NOS_*`  | tumor cells of PTCL-NOS                         | 9  |
| `EBV_*`  | tumor cells of EBV+ PTCL                        | 2  |

Reference group for every comparison: normal CD4 T cells pooled from 6 reactive lymph nodes
(sn04, sn08, sn09, sn10, sn15, sn16; 11,371 cells).

Processed data are not included. Only code is provided.

## Layout

```
01_define_groups/define_de_groups.py      tumor / CD4 cell groups -> de_groups.csv
02_differential_expression/run_mast.R     Seurat FindMarkers (MAST), one tumor group vs CD4
02_differential_expression/run_all.sh     all 29 comparisons
03_fgsea/run_fgsea.R                      fgsea on one DE table, MSigDB H and C1-C7
03_fgsea/run_all.sh                       all 29 DE tables
04_figure3a/samples.R                     sample names and order
04_figure3a/selected_pathways.txt         the 29 gene sets shown, in two blocks
04_figure3a/candidate_pathways.R          candidate gene sets the 29 were chosen from
04_figure3a/build_df_fig3a.R              collects fgsea results of the 29 gene sets -> df_Fig3A.csv
04_figure3a/plot_Fig3A.R                  df_Fig3A.csv -> plot_Fig3A_v3.pdf
```

## Requirements

* Python 3 with `pandas` (step 1 reads only metadata tables)
* R 4.4 with `Seurat` v5, `MAST`, `anndata` (R package, via reticulate), `fgsea`, `msigdbr`, `tidyverse`

Versions used for the checks below: R 4.4.3, Seurat 5.5.0, MAST 1.32.0, fgsea 1.32.4, ggplot2 4.0.0.
Seurat must be v5: FindMarkers was run with the v5 default `logfc.threshold = 0.1`, and the
original DE tables contain only genes with |avg_log2FC| > 0.1.

**MSigDB version.** The original fgsea runs used two msigdbr versions (see
[Note on sample sn13](#note-on-sample-sn13)):

* 28 samples (Mar 2025): **msigdbr 7.5.1** (MSigDB v7.5.1)
* `AITL_sn13` (re-run Aug 2025): **msigdbr 25.1.1** (MSigDB 2025.1; adds `KEGG_MEDICUS_*` sets and renames some GO terms)

To reproduce the original results, install msigdbr 7.5.1
(`install.packages('https://cran.r-project.org/src/contrib/Archive/msigdbr/msigdbr_7.5.1.tar.gz', repos = NULL)`).

The original jobs ran on SLURM (24 h, 50 GB RAM, 5 CPUs per job).

## Upstream inputs (not included)

`define_de_groups.py --data-dir <DIR>` and `run_mast.R` use these outputs of the upstream
snRNA-seq integration and sub-clustering analysis. They are the same inputs as the CellPhoneDB analysis:

| File | Content used |
|---|---|
| `Integration/raw_cellbender_DF_cat.h5ad` | CellBender-corrected, doublet-filtered raw counts of all samples. Cell barcodes are prefixed with the sample ID (`snXX_`). |
| `Integration/meta.csv`        | `leiden_res0.7`: clusters of the integrated object (major cell types) |
| `Subtype/T-cell/meta.csv`     | `leiden_res0.7`, `batch`: T-cell sub-clusters |
| `Subtype/leiden_14/meta.csv`  | `leiden_res0.9`: mast cells (10) and low-quality sub-clusters to remove (4, 5, 7, 11) |
| `Subtype/leiden_7/meta.csv`   | `leiden_res0.8`: low-quality sub-clusters to remove (2, 7, 10, 11) |
| `Subtype/Macrophage/meta.csv` | `leiden_res0.1`: cDC (1) |

## How to run

```bash
# 1. Cell groups (tumor cells of each sample; pooled normal CD4 T cells)
python 01_define_groups/define_de_groups.py --data-dir <DIR> --out de_groups.csv

# 2. Differential expression, tumor vs CD4 (29 comparisons; split into jobs if needed)
bash 02_differential_expression/run_all.sh <DIR>/Integration/raw_cellbender_DF_cat.h5ad de_groups.csv de_output

# 3. fgsea
bash 03_fgsea/run_all.sh de_output fgsea_output

# 4. Figure 3A
Rscript 04_figure3a/candidate_pathways.R fgsea_output candidate_pathways.csv   # optional
Rscript 04_figure3a/build_df_fig3a.R     fgsea_output df_Fig3A.csv
Rscript 04_figure3a/plot_Fig3A.R         df_Fig3A.csv plot_Fig3A_v3.pdf
```