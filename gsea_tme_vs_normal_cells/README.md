# TME cells of TFHL and PTCL-NOS vs normal lymph nodes: differential expression and GSEA (snRNA-seq)

This folder contains the code that produces the dot plot **`plot_fgsea_tme_condensed_v2_rd.pdf`**.
For four tumor-microenvironment (TME) cell types, the cells from lymphoma samples are compared with
the same cell type from reactive lymph nodes. The comparison is a differential expression test
(MAST) followed by pre-ranked GSEA (fgsea). The plot shows 16 selected gene sets. Each dot is coloured
by NES and sized by −log10(padj), and only gene sets with padj < 0.05 are drawn.

| Column | Lymphoma cells | Reference cells |
|---|---|---|
| `<CellType>-AITL` | that cell type pooled from 12 nodal TFH lymphoma (AITL / TFHL) samples | that cell type pooled from 6 reactive lymph nodes |
| `<CellType>-NOS`  | that cell type pooled from 9 PTCL-NOS samples                         | same |

Cell types: `B-cell`, `Plasma` (left panel), `cDC`, `Macrophage` (right panel). This gives 8 comparisons.
Reactive lymph nodes: sn04, sn08, sn09, sn10, sn15, sn16.

Processed data are not included. Only code is provided.

## Layout

```
01_define_groups/define_de_groups.py          TME cell groups (4 cell types x Normal/AITL/NOS) -> de_groups.csv
02_differential_expression/run_mast.R         Seurat FindMarkers (MAST), one cell type, lymphoma vs normal
02_differential_expression/run_all.sh         all 8 comparisons
03_fgsea/run_fgsea.R                          fgsea on one DE table, MSigDB H and C1-C7
03_fgsea/run_all.sh                           all 8 DE tables
04_figure/comparisons.R                       comparison names, order and panels
04_figure/selected_pathways.txt               the 16 gene sets shown, top to bottom
04_figure/candidate_pathways.R                candidate gene sets the 16 were chosen from
04_figure/build_df_fig.R                      collects fgsea results of the 16 -> df_fgsea_tme_condensed.csv
04_figure/plot_fgsea_tme_condensed.R          df_fgsea_tme_condensed.csv -> plot_fgsea_tme_condensed_v2_rd.pdf
```

## Requirements

* Python 3 with `pandas` (step 1 reads only metadata tables)
* R 4.4 with `Seurat` v5, `MAST`, `anndata` (R package, via reticulate), `fgsea`, `msigdbr`,
  `tidyverse`, `patchwork`, `cowplot`

Versions used for the checks below: R 4.4.3, Seurat 5.5.0, MAST 1.32.0, fgsea 1.32.4, ggplot2 4.0.0,
patchwork 1.3.2, cowplot 1.2.0.
Seurat must be v5: FindMarkers was run with the v5 default `logfc.threshold = 0.1`.

**MSigDB version.** All original fgsea runs (Mar 2025) used **msigdbr 7.5.1** (MSigDB v7.5.1; no
`KEGG_MEDICUS_*` sets). Newer msigdbr releases rename or drop some GO terms. To reproduce the original
results, install msigdbr 7.5.1
(`install.packages('https://cran.r-project.org/src/contrib/Archive/msigdbr/msigdbr_7.5.1.tar.gz', repos = NULL)`).

## Upstream inputs (not included)

These are the same inputs as the tumor-vs-CD4 GSEA and the CellPhoneDB analysis:

| File | Content used |
|---|---|
| `Integration/raw_cellbender_DF_cat.h5ad` | CellBender-corrected, doublet-filtered raw counts of all samples. Cell barcodes are prefixed with the sample ID (`snXX_`). |
| `Integration/meta.csv`        | `leiden_res0.7`: clusters of the integrated object (major cell types) |
| `Subtype/T-cell/meta.csv`     | `leiden_res0.7`: T-cell sub-clusters (they do not change the four TME cell types; kept so the rules match the original label files) |
| `Subtype/leiden_14/meta.csv`  | `leiden_res0.9`: mast cells (10) and low-quality sub-clusters to remove (4, 5, 7, 11) |
| `Subtype/leiden_7/meta.csv`   | `leiden_res0.8`: low-quality sub-clusters to remove (2, 7, 10, 11) |
| `Subtype/Macrophage/meta.csv` | `leiden_res0.1`: cDC (1) |

## How to run

```bash
# 1. Cell groups
python 01_define_groups/define_de_groups.py --data-dir <DIR> --out de_groups.csv

# 2. Differential expression, lymphoma vs normal (8 comparisons)
bash 02_differential_expression/run_all.sh <DIR>/Integration/raw_cellbender_DF_cat.h5ad de_groups.csv de_output

# 3. fgsea
bash 03_fgsea/run_all.sh de_output fgsea_output

# 4. Figure
Rscript 04_figure/candidate_pathways.R        fgsea_output candidate_pathways.csv   # optional
Rscript 04_figure/build_df_fig.R              fgsea_output df_fgsea_tme_condensed.csv
Rscript 04_figure/plot_fgsea_tme_condensed.R  df_fgsea_tme_condensed.csv plot_fgsea_tme_condensed_v2_rd.pdf
```