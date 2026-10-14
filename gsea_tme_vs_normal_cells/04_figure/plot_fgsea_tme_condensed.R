# Dot plot of fgsea results for 16 gene sets in TME cells of TFHL (AITL) and PTCL-NOS vs the same
# cell types in reactive lymph nodes. Concise version of ntbk_fgsea_tme_condensed.ipynb
# (the cells writing plot_fgsea_tme_condensed_v2_rd.pdf).
#
# Usage: Rscript plot_fgsea_tme_condensed.R [df_fgsea_tme_condensed.csv] [plot_fgsea_tme_condensed_v2_rd.pdf]

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(cowplot)
})

args <- commandArgs(TRUE)
IN_CSV  <- if (length(args) >= 1) args[1] else 'df_fgsea_tme_condensed.csv'
OUT_PDF <- if (length(args) >= 2) args[2] else 'plot_fgsea_tme_condensed_v2_rd.pdf'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

source(file.path(HERE, 'comparisons.R'))
selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')
displayed_pathway_order <- rev(selected$pathway)   # bottom to top

# Colour: NES capped at +/-2.5. Size: -log10(padj), padj floored at 1e-10.
# Only gene sets with padj < 0.05 are drawn ('_rd' version).
df <- read.csv(IN_CSV)
df$NES <- pmin(pmax(df$NES, -2.5), 2.5)
df$padj[which(df$padj < 1e-10)] <- 1e-10
df <- df %>% filter(padj < 0.05)

make_dotplot <- function(subset_df, samples, show_y_labels = TRUE, show_legend = TRUE) {
    ggplot(subset_df, aes(x = Sample, y = pathway, fill = NES, size = -log10(padj))) +
        # stroke fixed to the ggplot2 3.x default (ggplot2 4 derives it from the cowplot theme)
        geom_point(shape = 21, color = 'black', stroke = 0.5, show.legend = show_legend) +
        cowplot::theme_cowplot() +
        scale_fill_gradient2(midpoint = 0, limits = c(-2.5, 2.5),
                             low = '#2166ac', mid = '#f7f7f7', high = '#b2182b', space = 'Lab') +
        scale_x_discrete(limits = samples) +
        scale_y_discrete(limits = displayed_pathway_order) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              text = element_text(size = 14),
              axis.text.y = if (show_y_labels) element_text() else element_blank(),
              legend.position = if (show_legend) 'right' else 'none') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
        labs(title = '') +   # empty title, kept because it sets the top margin of the original
        if (show_legend) guides(fill = guide_colorbar(title = 'NES'),
                                size = guide_legend(title = '-log10(padj)'))
        else guides(fill = 'none', size = 'none')
}

# As in the original, each panel trains its own size scale; the legend is the right panel's.
p_left  <- df %>% filter(Sample %in% PANEL_LEFT) %>%
    make_dotplot(PANEL_LEFT, show_y_labels = TRUE, show_legend = FALSE)
p_right <- df %>% filter(Sample %in% PANEL_RIGHT) %>%
    make_dotplot(PANEL_RIGHT, show_y_labels = FALSE, show_legend = TRUE)

p <- p_left + p_right + plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'right')

pdf(OUT_PDF, width = 8, height = 6)
print(p)
invisible(dev.off())
