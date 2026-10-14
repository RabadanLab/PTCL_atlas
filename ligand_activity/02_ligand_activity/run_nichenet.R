# Sender-agnostic NicheNet ligand activity of one ecosystem gene set in one receiver cell type.
# Replaces {M1,M2}_to_{macrophage,endothelium,stromal}/ntbk_*_fixed_goi_rebak.ipynb (and Prepare.R)
# of the original analysis, which ran the same steps once per run.
#
# Usage: Rscript run_nichenet.R <raw_cellbender_DF_cat.h5ad> <meta_annotation_Lin.csv> \
#                               <genesets_of_interest.csv> <nichenet_dir> <out_dir> [run ...]
#   genesets_of_interest.csv: 01_genesets/genesets_of_interest.csv
#   nichenet_dir: contains lr_network_human_21122021.rds and ligand_target_matrix_nsga2r_final.rds
#                 (NicheNet v2 networks, https://zenodo.org/record/7074291)
#   run:          e.g. M1_to_macrophage (default: all 6)
# Output: <out_dir>/<run>.csv, the predict_ligand_activities table (test_ligand, auroc, aupr,
#         aupr_corrected, pearson), one row per potential ligand.

suppressPackageStartupMessages({
    library(nichenetr)
    library(Seurat)
    library(anndata)
    library(dplyr)
})

args <- commandArgs(TRUE)
H5AD     <- args[1]
META     <- args[2]
GENESETS <- args[3]
NN_DIR   <- args[4]
OUT_DIR  <- args[5]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

genesets <- read.csv(GENESETS)
to_run <- if (length(args) >= 6) args[-(1:5)] else unique(genesets$run)

lr_network <- readRDS(file.path(NN_DIR, 'lr_network_human_21122021.rds')) %>% distinct(from, to)
ligand_target_matrix <- readRDS(file.path(NN_DIR, 'ligand_target_matrix_nsga2r_final.rds'))

# Cell-type labels of all samples (the `reanno` column). The receiver is taken from all
# samples, not only those of one ecosystem.
meta <- read.csv(META, row.names = 1)
receivers <- unique(genesets$receiver[genesets$run %in% to_run])
cells <- rownames(meta)[meta$reanno %in% receivers]

# Raw (CellBender-corrected) counts. LogNormalize is per cell, so normalising after
# subsetting to the receiver cells gives the same values as normalising the full object.
adata <- read_h5ad(H5AD)
cells <- intersect(adata$obs_names, cells)
adata <- adata[cells, ]
mydata <- CreateSeuratObject(counts = t(adata$X))
rm(adata); invisible(gc())
mydata <- NormalizeData(mydata, normalization.method = 'LogNormalize', scale.factor = 10000)
Idents(mydata) <- meta[colnames(mydata), 'reanno']
print(table(Idents(mydata)))

for (RUN in to_run) {
    receiver <- genesets$receiver[genesets$run == RUN][1]
    geneset_oi <- genesets$gene[genesets$run == RUN]

    # Genes expressed in >= 5% of the receiver cells, their receptors, and all ligands of
    # those receptors (identical for the M1 and M2 runs of one receiver).
    expressed_genes_receiver <- get_expressed_genes(receiver, mydata, pct = 0.05)
    expressed_receptors <- intersect(unique(lr_network$to), expressed_genes_receiver)
    potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
    background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    cat(RUN, ': receiver', receiver, '| expressed genes', length(expressed_genes_receiver),
        '| receptors', length(expressed_receptors), '| potential ligands', length(potential_ligands),
        '| background', length(background_expressed_genes), '| gene set', length(geneset_oi), '\n')

    ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                   background_expressed_genes = background_expressed_genes,
                                                   ligand_target_matrix = ligand_target_matrix,
                                                   potential_ligands = potential_ligands)
    write.csv(ligand_activities, file.path(OUT_DIR, paste0(RUN, '.csv')))
}
