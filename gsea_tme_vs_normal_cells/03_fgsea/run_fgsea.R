# Pre-ranked GSEA (fgsea) of one TME-vs-normal comparison against MSigDB collections.
# Same as {Macrophage,Plasma,B-cell,cDC}/{AITL,NOS}_vs_Normal/runfGSEA.R of the original analysis.
#
# Usage: Rscript run_fgsea.R <de.csv> <out_dir>
#   <de.csv>: FindMarkers table from 02_differential_expression/run_mast.R
# Output: <out_dir>/set_<collection>.Rdata, each holding the fgsea result table `fgseaRes`.

suppressPackageStartupMessages({
    library(fgsea)
    library(msigdbr)
    library(tidyverse)
})

args <- commandArgs(TRUE)
DE_CSV  <- args[1]
OUT_DIR <- args[2]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df_de <- read.csv(DE_CSV, row.names = 1)
print(dim(df_de))

# Genes ranked by avg_log2FC (lymphoma vs normal); mitochondrial genes removed.
# Only genes reported by FindMarkers are ranked (min.pct 0.01, |avg_log2FC| > 0.1).
ranks <- df_de %>% arrange(desc(avg_log2FC))
ranks <- setNames(ranks$avg_log2FC, rownames(ranks))
ranks <- ranks[-grep('MT-', names(ranks))]

# The figure uses H, C2 (KEGG) and C5 (GO:BP); all collections were run originally.
for (LABEL in c('H', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7')) {
    msigdbr_df <- msigdbr(species = 'Homo sapiens', category = LABEL)
    fgsea_sets <- msigdbr_df %>% split(x = .$gene_symbol, f = .$gs_name)
    # Unlike the tumor-vs-CD4 runs, nPermSimple = 100000 and eps = 0 (exact small p-values).
    # fgseaMultilevel otherwise; no seed was set in the original runs.
    fgseaRes <- fgsea(fgsea_sets, stats = ranks, minSize = 15, maxSize = 500,
                      nPermSimple = 100000, eps = 0)
    save(fgseaRes, file = file.path(OUT_DIR, paste0('set_', LABEL, '.Rdata')))
}
