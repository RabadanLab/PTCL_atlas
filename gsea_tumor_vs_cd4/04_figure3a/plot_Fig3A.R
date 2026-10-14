# Fig. 3A heatmap: signed fgsea significance of 29 gene sets in tumor cells vs normal CD4 T cells,
# one column per tumor sample. Concise version of ntbk_Fig3A.ipynb (the cells writing plot_Fig3A_v3.pdf).
#
# Usage: Rscript plot_Fig3A.R [df_Fig3A.csv] [plot_Fig3A_v3.pdf]

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
IN_CSV  <- if (length(args) >= 1) args[1] else 'df_Fig3A.csv'
OUT_PDF <- if (length(args) >= 2) args[2] else 'plot_Fig3A_v3.pdf'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')
common_pathways    <- selected$pathway[selected$group == 'common']
aitl_high_pathways <- selected$pathway[selected$group == 'aitl_high']

df <- read.csv(IN_CSV)
df$Diagnosis <- sapply(strsplit(df$Sample, '_'), `[`, 1)

# significance = -log10(padj) signed by NES; padj floored at 1e-10.
# Gene sets whose padj is NA (fgsea could not estimate it) are shown in grey.
df$padj[which(df$padj < 1e-10)] <- 1e-10
df$significance <- (-log10(df$padj)) * sign(df$NES)

# pathway x sample matrices, used only for ordering (untested -> 0, padj NA -> NA)
mat_aitl <- xtabs(significance ~ pathway + Sample, data = df %>% filter(Diagnosis == 'AITL')) %>% as.matrix
mat_nos  <- xtabs(significance ~ pathway + Sample, data = df %>% filter(Diagnosis == 'NOS')) %>% as.matrix

# Columns: within TFHL (AITL) and within PTCL-NOS, the column order of hierarchical clustering
# (Euclidean distance, complete linkage; the pheatmap defaults used originally), reversed.
# EBV+ samples last.
col_order <- function(mat) rev(colnames(mat)[hclust(dist(t(mat)), method = 'complete')$order])
reordered_samples <- c(col_order(mat_aitl), col_order(mat_nos), 'EBV_sn38', 'EBV_sn47')

# Rows: two blocks (common, aitl_high), each sorted by mean significance over the
# TFHL + PTCL-NOS samples, highest at the top.
mat_all <- cbind(mat_aitl, mat_nos)
row_order <- function(pathways) {
    m <- rowMeans(mat_all[pathways, , drop = FALSE], na.rm = TRUE)
    rev(pathways[order(m)])
}
reordered_pathways <- c(row_order(common_pathways), row_order(aitl_high_pathways))

# Colour scale capped at +5. As in the original, only positive values are capped; the one value
# below -5 (KEGG_CELL_ADHESION_MOLECULES_CAMS in NOS_sn32, -5.03) is kept as is.
df$significance[df$significance > 5] <- 5

p <- ggplot(df, aes(x = Sample, y = pathway, fill = significance)) +
    geom_tile() +
    scale_y_discrete(limits = rev(reordered_pathways)) +
    scale_x_discrete(limits = reordered_samples) +
    scale_fill_gradient2(low = '#2166ac', mid = '#f7f7f7', high = '#b2182b', midpoint = 0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = 'black'))

pdf(OUT_PDF, width = 12, height = 7)
print(p)
invisible(dev.off())
