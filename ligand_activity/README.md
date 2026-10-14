# NicheNet ligand activity of TME programs in two lymphoma ecosystems (snRNA-seq)

This folder contains the code that produces the three scatter plots
**`plot_ligand_activity_nq_{macro,endothelium,stromal}_fixed_goi.pdf`**.
For three receiver cell types of the tumor microenvironment (macrophage, endothelium, stromal), NicheNet
(sender-agnostic ligand activity) asks which ligands best explain a gene program of that cell type.
Two programs are tested per cell type, one for ecosystem M1 (ecosystem 2) and one for ecosystem M2
(ecosystem 6). Each program is the top-gene list of a Spectra factor. Each plot compares the two runs of
one receiver. Every dot is a potential ligand at (activity for M1, activity for M2), after quantile
normalisation of the two runs. The top 10 ligands of each run are labelled in red.

| Plot | Receiver | M1 (ecosystem 2) gene set | M2 (ecosystem 6) gene set |
|---|---|---|---|
| `plot_ligand_activity_nq_macro_fixed_goi.pdf`       | macrophage  | factor `35-X-Mac-X-35` | factor `37-X-Mac-X-37` |
| `plot_ligand_activity_nq_endothelium_fixed_goi.pdf` | endothelium | factor `73-X-endo-X-73` | factor `74-X-endo-X-74` |
| `plot_ligand_activity_nq_stromal_fixed_goi.pdf`     | stromal     | factors `87-X-fibro-X-fibro_EGF_response` + `93-X-fibro-X-93` | factor `91-X-fibro-X-fibro_TGFB1_response` |

Processed data are not included. Only code is provided.

## Layout

```
01_genesets/runs.R                       the 6 runs: receiver and Spectra factor(s) of each
01_genesets/extract_genesets.R           Spectra top-gene table -> genesets_of_interest.csv
01_genesets/genesets_of_interest.csv     the 6 gene sets of interest (output of extract_genesets.R)
02_ligand_activity/run_nichenet.R        NicheNet ligand activity, one receiver x one gene set
02_ligand_activity/run_all.sh            all 6 runs
03_figure/plot_ligand_activity_nq.R      M1 and M2 runs of one receiver -> scatter plot + table
03_figure/run_all.sh                     all 3 receivers
```

## Requirements

* R 4.4 with `nichenetr` (v2), `Seurat` v5, `anndata` (R package, via reticulate), `tidyverse`,
  `ggrepel`, `preprocessCore`

Versions used for the checks below: R 4.4.3, nichenetr 2.2.1.1, Seurat 5.5.0, ggplot2 4.0.0,
ggrepel 0.9.6, preprocessCore 1.71.2.

## Inputs (not included)

| File | Content used |
|---|---|
| `Integration/raw_cellbender_DF_cat.h5ad` | CellBender-corrected, doublet-filtered raw counts of all samples. Cell barcodes are prefixed with the sample ID (`snXX_`). |
| `AnnotationLin/meta_annotation_Lin.csv`  | `reanno`: cell-type label of each cell (`macrophage`, `endothelium`, `stromal`, ...) |
| `topgenes_per_factor.csv`                | Spectra result: one column per factor, top 200 genes ranked by loading. Only needed to re-make `genesets_of_interest.csv`. |
| NicheNet v2 human networks ([Zenodo 7074291](https://zenodo.org/record/7074291)) | `lr_network_human_21122021.rds`, `ligand_target_matrix_nsga2r_final.rds` |

The NicheNet files can be downloaded with:

```bash
mkdir -p nichenet
for f in lr_network_human_21122021.rds ligand_target_matrix_nsga2r_final.rds; do
    wget -O nichenet/$f https://zenodo.org/record/7074291/files/$f
done
```

## How to run

```bash
# 1. Gene sets of interest (optional: the result is included)
Rscript 01_genesets/extract_genesets.R topgenes_per_factor.csv 01_genesets/genesets_of_interest.csv

# 2. NicheNet ligand activity (6 runs)
bash 02_ligand_activity/run_all.sh <DIR>/Integration/raw_cellbender_DF_cat.h5ad \
    <DIR>/AnnotationLin/meta_annotation_Lin.csv nichenet activity_output

# 3. Figures
bash 03_figure/run_all.sh activity_output figure_output
```

Step 2 takes about TIME_PLACEHOLDER and RAM_PLACEHOLDER of RAM.

## Method details

**Gene sets of interest.** For each run, the top 50 genes of the Spectra factor(s) in `runs.R`. For
M1_to_stromal, the top 50 genes of two factors are combined (89 unique genes). The other runs have 50 genes.
Genes that are not expressed in the receiver (see below) are outside the background and do not count.

**Receiver cells.** All cells of all samples whose `reanno` label is the receiver cell type:
macrophage 16,689, endothelium 2,869, stromal 8,262 cells. The ecosystem is not used to select cells.
The two runs of one receiver use the same cells and differ only in the gene set of interest.

**Ligand activity.** Raw counts are log-normalised (`NormalizeData`, LogNormalize, scale factor 10^4).
Then, as in the NicheNet v2 sender-agnostic vignette:

* expressed genes: `get_expressed_genes(receiver, pct = 0.05)`, i.e. detected in ≥ 5% of the receiver cells;
* potential ligands: all ligands in `lr_network` (distinct `from`–`to` pairs) whose receptor is an expressed gene;
* background: expressed genes that are in the ligand–target matrix;
* `predict_ligand_activities(geneset, background_expressed_genes, ligand_target_matrix, potential_ligands)`.

| Receiver | Expressed genes | Expressed receptors | Potential ligands | Background genes |
|---|---|---|---|---|
| macrophage  | 6,879 | 252 | 760 | 6,369 |
| endothelium | 7,122 | 297 | 822 | 6,587 |
| stromal     | 6,899 | 274 | 798 | 6,371 |

**Scatter plot.**

* Values: `aupr_corrected` (AUPR minus the AUPR of a random predictor) of each potential ligand in the M1
  and M2 runs of the receiver. The two vectors are quantile-normalised together
  (`preprocessCore::normalize.quantiles`) so that both axes have the same distribution.
* Red, labelled: the top 10 ligands of each run by `aupr_corrected` (before normalisation). For the
  endothelium, TGFB1 is labelled in addition. This gives 19 (macrophage), 21 (endothelium) and 20
  (stromal) labelled ligands.
* Axes: 0 to the maximum normalised activity + 0.005; dashed line y = x.
* Size 6 × 6 in.
* Label positions come from `ggrepel`. The notebook set no seed, so they differed between runs; the
  script sets `seed = 1`.

The table written with each plot (`aupr_corrected_<receiver>.csv`) has the `aupr_corrected` values
of both runs before (`aupr_corrected_ecosystem_{2,6}`) and after (`..._nq`) quantile normalisation.

## Correspondence to the original working files

Paths are relative to `snRNA/Nichenet/`.

| This repository | Original files |
|---|---|
| `runs.R`, `extract_genesets.R` | `geneset_oi` cell of `{M1,M2}_to_{macrophage,endothelium,stromal}/ntbk_*_fixed_goi_rebak.ipynb`; macrophage factors from `mac_spectra_top50_genes.csv`, the others from `topgenes_per_factor.csv` |
| `run_nichenet.R`, `run_all.sh` | `M1_to_macrophage/Prepare.R` → `obj_normalized_all_samples.rds`; `{M1,M2}_to_*/ntbk_*_fixed_goi_rebak.ipynb` (Nov 2025) → `ligand_activities_agnostic_fixed_goi.csv` |
| `plot_ligand_activity_nq.R`, `run_all.sh` | `Summary/ntbk_regulatory_potential_{macrophage,endothelial,stromal}_fixed_goi.ipynb` (last cells) → `plot_ligand_activity_nq_{macro,endothelium,stromal}_fixed_goi.pdf`, `aupr_corrected_{microphage,endothelial,stromal}.csv` |

`mac_spectra_top50_genes.csv` is identical to the first 50 rows of the macrophage factor columns of
`topgenes_per_factor.csv`, so one file is used here. The notebook that made `topgenes_per_factor.csv`
(Spectra) is not in this project folder.

Earlier versions, not used for these plots:

* `{M1,M2}_to_*/ntbk_{M1,M2}_to_*.ipynb` (Apr–May 2025): the gene set of interest was the MAST DE genes of
  the receiver in the ecosystem samples vs reactive lymph nodes (`geneset_oi.txt`), and the receiver cells came
  only from those samples (`obj_{M1,M2}_normal.rds`, made by `Prepare.R` from the sample lists
  M1: sn01, sn02, sn03, sn12, sn31, sn32, sn35, sn38, sn41, sn45; M2: sn33, sn36, sn39, sn44, sn47, sn49).
  → `ligand_activities_agnostic.csv`. The macrophage version of the scatter plot made from these runs is
  `Summary/plot_ligand_activity_nq.pdf` (`Summary/ntbk_macrophage_hubs.ipynb`).
* `{M1,M2}_to_*/ntbk_*_fixed_goi.ipynb` (Jun 2025): Spectra gene sets, but receiver cells still only
  from the ecosystem and normal samples. Its `ligand_activities_agnostic_fixed_goi.csv` was overwritten
  by the `_rebak` notebooks in Nov 2025.
* `Summary/ntbk_regulatory_potential_{macrophage,endothelial,stromal}.ipynb` (May–Jun 2025): ligand–target
  heatmaps only (`plot_ligand_target_mat_{macro,endo,stromal}_{M1,M2}.pdf`), no scatter plot.

Not included here: the ligand–target heatmaps `plot_ligand_target_mat_*_fixed_goi.pdf` made by the first
part of the same `Summary` notebooks (regulatory potential of the labelled ligands on the gene-set genes),
and the expression dot plots (`ntbk_check_expression_agnostic_*.ipynb`, `ntbk_check_ligand_expression_MRC1_RASSF4.ipynb`).

## Checks against the original intermediate files

CHECKS_PLACEHOLDER
