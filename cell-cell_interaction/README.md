# CellPhoneDB analysis: tumor T cell – microenvironment interactions (snRNA-seq)

This folder contains the code that produces the heatmap **`plot_heatmap_perc_sgnf.pdf`**.
For each selected ligand–receptor pair and each tumor-microenvironment (TME) cell type, the
heatmap shows the fraction of samples in which CellPhoneDB calls the interaction significant
(p < 0.05) between the TME cell type and one of these T-cell populations:

| Column | T-cell population | Samples (CellPhoneDB runs) |
|---|---|---|
| `CD4`  | normal CD4 T cells                        | 6 reactive lymph nodes pooled (sn04, sn08, sn09, sn10, sn15, sn16), 1 run |
| `TFH`  | normal T follicular helper cells          | same run as `CD4` |
| `TFHL` | tumor cells of nodal TFH lymphoma (AITL)  | 18 samples, 1 run per sample |
| `NOS`  | tumor cells of PTCL-NOS                   | 9 samples, 1 run per sample |

Processed data are not included. Only code is provided.

## Layout

```
01_prepare_inputs/prepare_cpdb_inputs.py   normalised counts + cell-type labels for CellPhoneDB
02_run_cellphonedb/run_cellphonedb.py      runs one CellPhoneDB statistical analysis
02_run_cellphonedb/run_all.sh              runs all 28 analyses (1 normal + 18 TFHL + 9 NOS)
03_heatmap/selected_pairs.txt              the 49 ligand-receptor pairs shown
03_heatmap/selected_pairs_mirror.txt       the same pairs with partners swapped (line-matched)
03_heatmap/build_df_plot.R                 collects p-values/means of the selected pairs -> df_plot.csv
03_heatmap/plot_heatmap_perc_sgnf.R        df_plot.csv -> plot_heatmap_perc_sgnf.pdf
```

## Requirements

* Python 3 with `scanpy`, `anndata`, `pandas`, `numpy`
* [CellPhoneDB](https://github.com/ventolab/CellphoneDB) v5 (the code was checked against v5.0.1)
  with database **cellphonedb-data v4.1.0** (`cellphonedb.zip`)
* R (≥ 4.0) with `tidyverse` (`dplyr`, `ggplot2`)

The original CellPhoneDB jobs ran on SLURM with 48 h, 100 GB RAM and 10 CPUs per job.

## Upstream inputs (not included)

`prepare_cpdb_inputs.py --data-dir <DIR>` expects these outputs of the upstream
snRNA-seq integration and sub-clustering analysis:

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
# 1. Inputs: normalize_total(1e4) + log1p, cell-type labels, split per sample
python 01_prepare_inputs/prepare_cpdb_inputs.py --data-dir <DIR> --out-dir cpdb_input

# 2. CellPhoneDB (28 runs; submit them as separate cluster jobs if needed)
bash 02_run_cellphonedb/run_all.sh cpdb_input cpdb_output path/to/cellphonedb.zip

# 3. Collect the selected pairs and plot
Rscript 03_heatmap/build_df_plot.R cpdb_output df_plot.csv
Rscript 03_heatmap/plot_heatmap_perc_sgnf.R df_plot.csv plot_heatmap_perc_sgnf.pdf
```

The output folder layout used by step 3 is:
```
cpdb_output/normal/{pvalues,means}.txt
cpdb_output/AITL/<sample>/{pvalues,means}.txt
cpdb_output/NOS/<sample>/{pvalues,means}.txt
```