# Differential expression of each tumor sample vs pooled normal CD4 T cells (Seurat FindMarkers, MAST).
# Replaces AITL_vs_CD4/RunMAST.R, NOS_vs_CD4/RunMAST.R and NOS_vs_CD4/RunMASTebv.R of the original
# analysis, which ran the same test once per tumor sample.
#
# Usage: Rscript run_mast.R <raw_cellbender_DF_cat.h5ad> <de_groups.csv> <out_dir> [group ...]
#   de_groups.csv: output of 01_define_groups/define_de_groups.py
#   group:         tumor groups to test (default: all groups other than CD4)
# Output: <out_dir>/<group>.csv, the FindMarkers table (tumor = ident.1, so avg_log2FC > 0 means
#         higher in tumor cells).

suppressPackageStartupMessages({
    library(Seurat)
    library(MAST)
    library(anndata)
})

args <- commandArgs(TRUE)
H5AD    <- args[1]
GROUPS  <- args[2]
OUT_DIR <- args[3]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

groups <- read.csv(GROUPS)
to_test <- if (length(args) >= 4) args[-(1:3)] else setdiff(unique(groups$group), 'CD4')

# Raw (CellBender-corrected) counts; LogNormalize is per cell, so normalising after
# subsetting to the cells used gives the same values as normalising the full object.
adata <- read_h5ad(H5AD)
adata <- adata[groups$cell, ]
mydata <- CreateSeuratObject(counts = t(adata$X))
rm(adata); invisible(gc())
mydata <- NormalizeData(mydata, normalization.method = 'LogNormalize', scale.factor = 10000)

cd4_cells <- groups$cell[groups$group == 'CD4']
for (GROUP in to_test) {
    # tumor cells first, then CD4 cells, as in the original scripts
    kept_cells <- c(groups$cell[groups$group == GROUP], cd4_cells)
    obj <- mydata[, kept_cells]
    Idents(obj) <- groups$group[match(colnames(obj), groups$cell)]
    print(table(Idents(obj)))

    # Defaults otherwise (logfc.threshold = 0.1 in Seurat v5, random.seed = 1 for the
    # down-sampling to 1000 cells per group).
    x <- FindMarkers(obj, ident.1 = GROUP, test.use = 'MAST', min.pct = 0.01,
                     max.cells.per.ident = 1000)
    write.csv(x, file.path(OUT_DIR, paste0(GROUP, '.csv')))
}
