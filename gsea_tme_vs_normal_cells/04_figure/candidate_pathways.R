# List the candidate gene sets from which the figure's pathways were chosen by hand.
# Concise version of the pathway-selection part of ntbk_fgsea_tme.ipynb
# (which wrote rslt_tme/plot_fgsea_tme_{H,C2,C5}.csv, the inputs of ntbk_fgsea_tme_condensed.ipynb).
#
# For each comparison and each of three gene-set groups (Hallmark; all of C2; GO:BP from C5),
# the 10 gene sets with the highest and the 10 with the lowest NES among those with padj < 0.05
# are taken; the candidates are the union over comparisons.
#
# Usage: Rscript candidate_pathways.R <fgsea_output_dir> [out.csv]

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
FGSEA_OUT <- if (length(args) >= 1) args[1] else 'fgsea_output'
OUT_CSV   <- if (length(args) >= 2) args[2] else 'candidate_pathways.csv'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

source(file.path(HERE, 'comparisons.R'))
selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')

# C2 is not restricted to KEGG here (the KEGG-only list, plot_fgsea_tme_kegg.csv, was not used)
gene_set_groups <- list(H = c('H', ''), C2 = c('C2', ''), GOBP = c('C5', 'GOBP'))

candidates <- NULL
for (COMPARISON in COMPARISONS) {
    for (GROUP in names(gene_set_groups)) {
        CATEGORY <- gene_set_groups[[GROUP]][1]
        PREFIX   <- gene_set_groups[[GROUP]][2]
        load(file.path(FGSEA_OUT, COMPARISON, paste0('set_', CATEGORY, '.Rdata')))
        sgnf <- fgseaRes %>% filter(padj < 0.05 & startsWith(pathway, PREFIX))
        # top_n keeps ties, as in the original notebook
        top <- c(top_n(sgnf, n = 10, wt = NES)$pathway, top_n(sgnf, n = -10, wt = NES)$pathway)
        if (length(top) > 0)
            candidates <- bind_rows(candidates, tibble(gene_set_group = GROUP, pathway = top, Sample = COMPARISON))
    }
}

candidates <- candidates %>%
    group_by(gene_set_group, pathway) %>%
    summarise(n_comparisons_top = n_distinct(Sample), .groups = 'drop') %>%
    mutate(in_figure = pathway %in% selected$pathway) %>%
    arrange(gene_set_group, desc(n_comparisons_top))

print(table(candidates$gene_set_group))
stopifnot(all(selected$pathway %in% candidates$pathway))
write.csv(candidates, OUT_CSV, row.names = FALSE)
