# Differential expression of each TME cell type in lymphoma (TFHL/AITL or PTCL-NOS) vs the same cell
# type in reactive lymph nodes (Seurat FindMarkers, MAST).
# Replaces {Macrophage,Plasma,B-cell,cDC}/{AITL,NOS}_vs_Normal/RunMAST.R of the original analysis,
# which ran the same test once per cell type and diagnosis.
#
# Usage: Rscript run_mast.R <raw_cellbender_DF_cat.h5ad> <de_groups.csv> <out_dir> [comparison ...]
#   de_groups.csv: output of 01_define_groups/define_de_groups.py
#   comparison:    <CellType>-<AITL|NOS>, e.g. Macrophage-AITL (default: all 8)
# Output: <out_dir>/<comparison>.csv, the FindMarkers table (lymphoma = ident.1, so avg_log2FC > 0
#         means higher in the lymphoma cells).

suppressPackageStartupMessages({
    library(Seurat)
    library(MAST)
    library(anndata)
})

# label in de_groups.csv -> name used for outputs and in the figure
CELL_TYPES <- c(macrophage = 'Macrophage', plasma = 'Plasma', `B-cell` = 'B-cell', cDC = 'cDC')

args <- commandArgs(TRUE)
H5AD    <- args[1]
GROUPS  <- args[2]
OUT_DIR <- args[3]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

groups <- read.csv(GROUPS, check.names = FALSE)
groups$cell_type <- CELL_TYPES[groups$cell_type]
all_comparisons <- as.vector(outer(CELL_TYPES, c('AITL', 'NOS'), paste, sep = '-'))
to_test <- if (length(args) >= 4) args[-(1:3)] else all_comparisons

# Raw (CellBender-corrected) counts; LogNormalize is per cell, so normalising after
# subsetting to the cells used gives the same values as normalising the full object.
adata <- read_h5ad(H5AD)
adata <- adata[groups$cell, ]
mydata <- CreateSeuratObject(counts = t(adata$X))
rm(adata); invisible(gc())
mydata <- NormalizeData(mydata, normalization.method = 'LogNormalize', scale.factor = 10000)

for (COMPARISON in to_test) {
    CELL_TYPE <- sub('-(AITL|NOS)$', '', COMPARISON)
    CONDITION <- sub('.*-', '', COMPARISON)
    # lymphoma cells first, then normal cells, as in the original scripts
    is_type <- groups$cell_type == CELL_TYPE
    kept_cells <- c(groups$cell[is_type & groups$condition == CONDITION],
                    groups$cell[is_type & groups$condition == 'Normal'])
    obj <- mydata[, kept_cells]
    Idents(obj) <- groups$condition[match(colnames(obj), groups$cell)]
    print(table(Idents(obj)))

    # Defaults otherwise (logfc.threshold = 0.1 in Seurat v5, random.seed = 1 for the
    # down-sampling to 1000 cells per group).
    x <- FindMarkers(obj, ident.1 = CONDITION, test.use = 'MAST', min.pct = 0.01,
                     max.cells.per.ident = 1000)
    write.csv(x, file.path(OUT_DIR, paste0(COMPARISON, '.csv')))
}
