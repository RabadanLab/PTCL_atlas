# TFHL tumor cells vs normal TFH cells: recurrently up-regulated genes (snRNA-seq)

This folder contains the code that produces the dot plot **`dotplot_recurrent_genes_TFHL_vs_TFH.pdf`**.
For each of 18 nodal TFH lymphoma (TFHL / AITL) samples, the tumor cells are compared with normal
TFH cells from reactive lymph nodes by differential expression (MAST). Genes significantly
up-regulated in many samples are then counted. The dot plot shows 23 recurrent genes (rows, grouped by
the number of samples in which they are significant) in normal TFH cells and TFHL tumor cells of
each sample (columns). Dot colour is mean expression scaled to 0–1 per gene, and dot size is the
fraction of cells expressing the gene.

| Group | Cells | Samples |
|---|---|---|
| `AITL_snXX` | tumor cells of one TFHL sample | 18 (sn01, sn02, sn03, sn05, sn06, sn07, sn12, sn13, sn14, sn29, sn31, sn35, sn36, sn39, sn43, sn45, sn46, sn49) |
| `TFH`       | normal TFH cells pooled from 6 reactive lymph nodes | sn04, sn08, sn09, sn10, sn15, sn16 (2,184 cells) |

Processed data are not included. Only code is provided.

## Layout

```
01_define_groups/define_de_groups.py           TFHL tumor / normal TFH cell groups -> de_groups.csv
02_differential_expression/run_mast.R          Seurat FindMarkers (MAST), one tumor sample vs TFH
02_differential_expression/run_all.sh          all 18 comparisons
03_recurrent_genes/recurrent_deg.R             significant genes per sample, count of samples per gene
04_dotplot/selected_genes.txt                  the 23 genes shown, top to bottom
04_dotplot/plot_dotplot_recurrent_genes.py     -> dotplot_recurrent_genes_TFHL_vs_TFH.pdf
```

## Requirements

* Python 3 with `pandas` (step 1), and `scanpy`, `anndata` (step 4)
* R 4.4 with `Seurat` v5, `MAST`, `anndata` (R package, via reticulate), `tidyverse`

Versions used for the checks below: R 4.4.3, Seurat 5.5.0, MAST 1.32.0; Python 3.12, scanpy 1.11.4,
anndata 0.12.2, matplotlib 3.10.6.
Seurat must be v5: FindMarkers was run with the v5 default `logfc.threshold = 0.1`, and the
original DE tables contain only genes with |avg_log2FC| > 0.1.

The original MAST jobs ran on SLURM (24 h, 50 GB RAM, 5 CPUs).

## Upstream inputs (not included)

`define_de_groups.py`, `run_mast.R` and `plot_dotplot_recurrent_genes.py` use these outputs of the
upstream snRNA-seq integration and sub-clustering analysis. They are the same inputs as the
tumor-vs-CD4 GSEA and the CellPhoneDB analysis:

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
# 1. Cell groups (tumor cells of each TFHL sample; pooled normal TFH cells)
python 01_define_groups/define_de_groups.py --data-dir <DIR> --out de_groups.csv

# 2. Differential expression, tumor vs TFH (18 comparisons; split into jobs if needed)
bash 02_differential_expression/run_all.sh <DIR>/Integration/raw_cellbender_DF_cat.h5ad de_groups.csv de_output

# 3. Recurrent genes -> deg_sgnf_TFHL_vs_TFH.csv, recurrent_deg_TFHL_vs_TFH.csv
Rscript 03_recurrent_genes/recurrent_deg.R de_output .

# 4. Dot plot
python 04_dotplot/plot_dotplot_recurrent_genes.py --data-dir <DIR> --out dotplot_recurrent_genes_TFHL_vs_TFH.pdf
```