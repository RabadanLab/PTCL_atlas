# Scatter plot of quantile-normalised NicheNet ligand activity, M1 (ecosystem 2) vs M2 (ecosystem 6),
# for one receiver cell type.
# Replaces the ligand-activity cells (last part) of Summary/ntbk_regulatory_potential_{macrophage,
# endothelial,stromal}_fixed_goi.ipynb.
#
# Usage: Rscript plot_ligand_activity_nq.R <activity_dir> <receiver> <out.pdf> <out.csv>
#   activity_dir: output of 02_ligand_activity (M1_to_<receiver>.csv, M2_to_<receiver>.csv)
#   receiver:     macrophage, endothelium or stromal
# Output: out.pdf, the scatter plot; out.csv, aupr_corrected of both runs before and after
#         quantile normalisation, one row per ligand.

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
    library(preprocessCore)
})

# Ligands labelled in addition to the top 10 of each run
EXTRA_LIGANDS <- list(macrophage = character(0), endothelium = 'TGFB1', stromal = character(0))
N_TOP <- 10

args <- commandArgs(TRUE)
ACT_DIR  <- args[1]
RECEIVER <- args[2]
OUT_PDF  <- args[3]
OUT_CSV  <- args[4]

ligand_activity_M1 <- read.csv(file.path(ACT_DIR, paste0('M1_to_', RECEIVER, '.csv')), row.names = 1)
ligand_activity_M2 <- read.csv(file.path(ACT_DIR, paste0('M2_to_', RECEIVER, '.csv')), row.names = 1)

# Labelled ligands: top 10 by aupr_corrected (before normalisation) of each run
M1_ligands <- ligand_activity_M1 %>% arrange(desc(aupr_corrected)) %>% pull(test_ligand) %>% head(N_TOP)
M2_ligands <- ligand_activity_M2 %>% arrange(desc(aupr_corrected)) %>% pull(test_ligand) %>% head(N_TOP)
ligands_to_show <- c(M1_ligands, M2_ligands, EXTRA_LIGANDS[[RECEIVER]]) %>% unique

# Both runs have the same receiver and hence the same potential ligands; the intersection is
# kept as in the notebook.
rownames(ligand_activity_M1) <- ligand_activity_M1$test_ligand
rownames(ligand_activity_M2) <- ligand_activity_M2$test_ligand
common <- intersect(rownames(ligand_activity_M1), rownames(ligand_activity_M2))
ligand_activity_M1 <- ligand_activity_M1[common, ]
ligand_activity_M2 <- ligand_activity_M2[common, ]

mat_nq <- cbind(ligand_activity_M1$aupr_corrected, ligand_activity_M2$aupr_corrected) %>%
    as.matrix %>% normalize.quantiles

df_activity <- data.frame(
    ligand = common,
    activity_nq_M1 = mat_nq[, 1],
    activity_nq_M2 = mat_nq[, 2]
)
df_activity$Highlight <- ifelse(df_activity$ligand %in% ligands_to_show, 'Yes', 'No')

MAX <- c(df_activity$activity_nq_M1, df_activity$activity_nq_M2) %>% max + 0.005

p <- ggplot(df_activity, aes(x = activity_nq_M1, y = activity_nq_M2)) +
    geom_point(aes(color = Highlight), size = 2, alpha = 0.5) +
    scale_color_manual(values = c('grey', 'red')) +
    theme_bw(base_size = 20) +
    geom_text_repel(
        data = subset(df_activity, Highlight == 'Yes'),
        aes(label = ligand),
        parse = TRUE,
        nudge_y = 0.002,
        size = 4.5,
        segment.size = 0.5,
        segment.color = 'black',
        force = 10,
        min.segment.length = 0,
        max.overlaps = Inf,
        # the notebook set no seed, so its label positions cannot be reproduced exactly
        seed = 1
    ) +
    xlim(c(0, MAX)) + ylim(c(0, MAX)) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', linewidth = 1) +
    theme(legend.position = 'none') +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black')
    )

pdf(OUT_PDF, width = 6, height = 6)
print(p)
invisible(dev.off())

df_out <- data.frame(
    aupr_corrected_ecosystem_2 = ligand_activity_M1$aupr_corrected,
    aupr_corrected_ecosystem_6 = ligand_activity_M2$aupr_corrected,
    aupr_corrected_ecosystem_2_nq = mat_nq[, 1],
    aupr_corrected_ecosystem_6_nq = mat_nq[, 2]
)
rownames(df_out) <- common
write.csv(df_out, OUT_CSV)
