# Heatmap of the fraction of samples in which each selected ligand-receptor pair is
# significant (CellPhoneDB p < 0.05) between the T-cell population of interest and each
# TME cell type, for normal CD4, normal TFH, TFHL (AITL) tumor cells and PTCL-NOS tumor cells.
# Concise version of ntbk_polished_dotplot_per_aitl.ipynb (heatmap cell) / ntbk_plot_heatmap.ipynb.
#
# Usage: Rscript plot_heatmap_perc_sgnf.R [df_plot.csv] [out.pdf] [--refined]
#   --refined: drop VCAM1_integrin_a4b1_complex and PDCD1_PDCD1LG2 and use the curated row
#              order from ntbk_plot_heatmap.ipynb (default reproduces plot_heatmap_perc_sgnf.pdf)

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
REFINED <- '--refined' %in% args
args <- setdiff(args, '--refined')
IN_CSV  <- if (length(args) >= 1) args[1] else 'df_plot.csv'
OUT_PDF <- if (length(args) >= 2) args[2] else 'plot_heatmap_perc_sgnf.pdf'

df_plot <- read.csv(IN_CSV, row.names = 1)

# Fraction of samples (unique celltype labels per diagnosis) with p < 0.05
n_samples <- df_plot %>% group_by(diagnosis) %>% summarise(N = n_distinct(celltype))
tmp <- df_plot %>%
    filter(diagnosis %in% c('CD4', 'TFH', 'TFHL', 'NOS')) %>%
    group_by(diagnosis, interacting_pair, tme) %>%
    summarise(num_sgnf = sum(p < 0.05), .groups = 'drop') %>%
    left_join(n_samples, by = 'diagnosis') %>%
    mutate(perc_sgnf = num_sgnf / N,
           tme = factor(tme, levels = c('CD4', 'CD8', 'Treg', 'macrophage', 'monocyte', 'cDC',
                                        'pDC', 'B.cell', 'GCB', 'plasma', 'stromal', 'endothelium')))

y_scale <- NULL
if (REFINED) {
    tmp <- tmp %>% filter(!interacting_pair %in% c('VCAM1_integrin_a4b1_complex', 'PDCD1_PDCD1LG2'))
    row_order <- c(
        'APP_CD74', 'CADM1_CADM1', 'CD28_CD80', 'CTLA4_CD80', 'TGFB1_TGFbeta_receptor1',
        'ANGPT2_TEK', 'ADGRE5_CD55', 'Dehydroepiandrosterone_bySTS_ESR1',
        'Dehydroepiandrosterone_bySTS_PPARA', 'CADM1_NECTIN3', 'ITGAL_ICAM1', 'CD226_NECTIN2',
        'TNFRSF10A_TNFSF10', 'CD46_JAG1', 'PDCD1_CD274', 'SIRPG_CD47', 'FGFR1_FGF2',
        'PDGFD_PDGFR_complex', 'TGFB1_TGFBR3', 'CLSTN1_NRXN3', 'CD40LG_CD40', 'SORL1_APP',
        'TIGIT_NECTIN2', 'SPN_ICAM1', 'SEMA4D_PLXNB2', 'CD44_TYROBP', 'SEMA4D_CD72',
        'ADGRE5_THY1', 'Glutamate_byGLS_and_SLC1A1_Glutamate_Kainate_1_4_complex',
        'ProstaglandinE2_byPTGES2_PTGER4', 'TNFRSF14_BTLA', 'NOTCH1_JAG1', 'VSIR_HLA-E',
        'CD38_PECAM1', 'TENM1_ADGRL1', 'SELL_SELPLG', 'APP_SORL1', 'CD28_CD86', 'CTLA4_CD86',
        'CD2_CD58', 'CD6_ALCAM', 'IGF1R_IGF1', 'TGFB1_TGFbeta_receptor2', 'TGFBR3_TGFB1',
        'TGFbeta_receptor1_TGFB1', 'TGFbeta_receptor2_TGFB1', 'ADGRE5_integrin_aVb3_complex')
    y_scale <- scale_y_discrete(limits = rev(row_order))
}

p_heatmap <- ggplot(tmp, aes(x = diagnosis, y = interacting_pair, fill = perc_sgnf)) +
    geom_tile() +
    facet_wrap(. ~ tme, nrow = 1) +
    scale_x_discrete(limits = c('CD4', 'TFH', 'TFHL', 'NOS')) +
    y_scale +
    scale_fill_gradient(low = 'gray', high = 'dark red') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab('') + ylab('')

pdf(OUT_PDF, width = 18, height = 10)
print(p_heatmap)
invisible(dev.off())
