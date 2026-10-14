# PTCL_atlas

This repository contains the custom code used to process, integrate and cluster the single-nucleus RNA-seq (snRNA-seq) data from 
the PTCL atlas. It also includes the down-sampling stability procedure used to choose the clustering resolution for each cell lineage.

## Contents

| Folder | Script | What it does | Manuscript section |
|---|---|---|---|
| `snRNA-seq_initial_clustering/` | `run_prepare.py` | Loads one sample's CellBender-filtered count matrix, keeps cells that passed DoubletFinder, and prefixes barcodes with the sample ID | Methods: *Preprocessing, integration, and clustering* |
| | `run_scanorama.py` | Normalises, selects genes highly variable in ≥10 samples, and integrates samples with Scanorama (`dimred = 50`) | Supplementary Methods: *Integration, batch effect removal and clustering* |
| | `run_clustering.py` | Builds the neighbour graph (30 Scanorama dimensions, 50 neighbours), UMAP, and Leiden clustering at 16 resolutions (0.05–1.5), with silhouette scores | Same |
| `snRNA-seq_subclustering/` | `SubClustering.py` | Re-integrates and re-clusters one lineage on the full data or on a random 90% down-sample | Supplementary Methods: stability-based sub-clustering |
| | `JaccardSimilarity.R` | Compares each down-sampled clustering with the full-data clustering (maximal Jaccard index) and computes the weighted stability score *S_R* across resolutions | Same (equations for *J* and *S_R*) |
| `Cell-cell_interaction/` | `ntbk_plot_heatmap.ipynb` | Summarises and plots CellPhoneDB (v3, statistical mode) ligand–receptor results for TFHL and PTCL, NOS versus normal tissue | Methods: *TME ecosystem and cell-cell interaction analysis* |
| `demo/` | `make_demo_data.py` | Generates a small **simulated** dataset (no patient data) to test the pipeline | — |

## System requirements

**Operating systems.**
Codes were run on an OS of Red Hat Enterprise Linux 9.6 with Python 3.12.11 and R 4.4.3.

## Demo

The patient data are under controlled access (see *Data availability* below), so the demo uses a **simulated** dataset that has the same structure as the real one. 
It has 12 samples × 300 nuclei × 3,000 genes, five cell types with distinct marker genes, and a sample-specific batch effect. 
At least 10 samples are needed because the pipeline keeps genes that are highly variable in ≥10 samples.

```bash
mkdir demo_run && cd demo_run && mkdir output
export MPLBACKEND=Agg   # run plots without a display

# 1. Simulate data (writes raw_cellbender_DF_cat.h5ad, raw_data.h5ad, random_select.pkl)
python ../demo/make_demo_data.py

# 2. Integration and initial clustering
python ../snRNA-seq_initial_clustering/run_scanorama.py     # -> integrated.h5ad
python ../snRNA-seq_initial_clustering/run_clustering.py    # -> integrated_processed.h5ad

# 3. Stability analysis for the sub-clustering resolution
python ../snRNA-seq_subclustering/SubClustering.py 1.0 1    # reference clustering on all cells
for i in $(seq 1 10); do
    python ../snRNA-seq_subclustering/SubClustering.py 0.9 $i   # ten 90% down-samples
done
Rscript ../snRNA-seq_subclustering/JaccardSimilarity.R      # -> plot_jaccard.pdf, plot_perc.pdf
```

**Expected output**

- `integrated.h5ad`: 3,600 nuclei with a 50-dimensional Scanorama embedding in `obsm["Scanorama"]`.
- `integrated_processed.h5ad`: the same object with a UMAP and Leiden clusters at 16 resolutions (`obs["leiden_res0.05"]` … `obs["leiden_res1.5"]`). On the demo, every resolution recovers the **5 simulated cell types** exactly. You can check this against `obs["true_celltype"]`.
- `output/subtypes_1.0_1.csv` and `output/subtypes_0.9_1.csv` … `subtypes_0.9_10.csv`: cluster labels at each resolution for the full data and each down-sample.
- `plot_jaccard.pdf`: the maximal Jaccard index of every cluster across down-samples, one page per resolution.
- `plot_perc.pdf`: the weighted stability score *S_R* versus resolution at four Jaccard thresholds (0.5, 0.6, 0.75, 0.9). Because the simulated cell types are well separated, *S_R* = 1.0 at every resolution on the demo. On real data *S_R* falls as resolution increases, and the optimal resolution is chosen where it starts to drop (Supplementary Methods).

**Expected run time on a normal desktop (2 cores):** about 10 minutes in total. That breaks down as about 5 s for simulation and integration, about 1 min for initial clustering, about 45 s for each `SubClustering.py` run (about 8 min for 11 runs), and about 5 s for the R script.

## Instructions for use

**Running on your own data**

1. Process each sample with Cell Ranger, then CellBender (default parameters), then DoubletFinder. For each sample, save the filtered CellBender matrix (`*_cellbender_output_filtered.h5`) and a one-column text file of barcodes that pass doublet removal (`kept_cells_DF.txt`).
2. Run `run_prepare.py <sample_folder> <h5_file>` once per sample. It writes `raw_cellbender_DF/<sample>.h5ad`. The script assumes the first four characters of the folder name are the sample ID (e.g. `sn01`).
   *Note:* the version used in the paper copies inputs from the lab's private AWS S3 bucket. If you are working with local files, remove the two `aws s3 cp` steps and place the files in the working directory.
3. Concatenate the per-sample files into `raw_cellbender_DF_cat.h5ad`, for example with `anndata.concat([...])`. Then apply per-sample QC filters for total UMI counts and mitochondrial content. The thresholds used for each sample are listed in Supplementary Table 2.
4. Run `run_scanorama.py`, then `run_clustering.py`, as in the demo. Choose the initial resolution from the silhouette plot. The paper used resolution 0.7.
5. For each lineage you want to sub-cluster, save its raw counts as `raw_data.h5ad` with a `batch` column. Create `random_select.pkl` with the cell masks: `demo/make_demo_data.py` contains a `make_random_select()` function that builds this file. Then run `SubClustering.py 1.0 1` and `SubClustering.py 0.9 <i>` for *i* = 1…100, followed by `JaccardSimilarity.R`.

## Specific analysis tasks

Please refer to the subfolders for codes and mannuals of specific analysis tasks, including differential gene expression, pathway enrichment, ligand activity, cell-cell interaction, etc.

## Data availability

snRNA-seq, bulk RNA-seq and targeted sequencing data are deposited in dbGaP under accession **phs004552.v1.p1** (controlled access).

## License

This code is released under the MIT License (see `LICENSE`).

## Contact

Teresa Palomero (tp2151@cumc.columbia.edu) · Raul Rabadan (rr2579@cumc.columbia.edu)
For questions about the code, please open an issue on this repository.

## Versions of tools and packages

| Tool | Version used in the paper |
| Cell Ranger | 7.0.0 |
| CellBender | 0.3.2 |
| DoubletFinder | 2.0.4 |
| Scanpy 1.11.4 | Scanorama 1.7.4 |
| Trimmomatic 0.39 | STAR 2.7.11b |
| HTSeq | 2.0.5 |
| TRUST4 | 1.1.6 |
| PolyPhen-2 | 2.2.3 |
| SIFT | 6.2.1 |
| REVEL | 1.3 |
| CNVkit | 0.9.12 |
| Numbat | 1.4.2 |
| MAST | 1.32.0 |
| Seurat | 5.3.0 |
| msigdbr | 25.1.1 |
| fgsea | 1.32.4 |
| MAGIC | 3.0.0 |
| UCell | 2.12.0 |
| Spectra | 1.0.1 |
| Cytopus | 0.3.0 |
| DESeq2 | 1.46.0 |
| CellPhoneDB | v3 |
| NicheNet | 2.2.0 |
| ssGSEA | 2.0 |
| GSEA/MSigDB | 25.1.1 |
| Prism | 10.5.0 |
| FlowJo | 10.9 |
