# List the candidate gene sets from which the Fig. 3A pathways were chosen by hand.
# Concise version of the pathway-selection part of ntbk_fgsea_tumor_vs_CD4.ipynb
# (which wrote rslt_tumor_vs_CD4/plot_fgsea_tumor_vs_CD4_{H,C5,kegg}.csv).
#
# For each sample and each of three gene-set groups (Hallmark; GO:BP from C5; KEGG from C2),
# the 10 gene sets with the highest and the 10 with the lowest NES among those with padj < 0.05
# are taken; the candidates are the union over samples.
#
# Usage: Rscript candidate_pathways.R <fgsea_output_dir> [out.csv]

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
FGSEA_OUT <- if (length(args) >= 1) args[1] else 'fgsea_output'
OUT_CSV   <- if (length(args) >= 2) args[2] else 'candidate_pathways.csv'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

source(file.path(HERE, 'samples.R'))
selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')

gene_set_groups <- list(H = c('H', ''), GOBP = c('C5', 'GOBP'), KEGG = c('C2', 'KEGG'))

candidates <- NULL
for (SAMPLE in SAMPLES) {
    for (GROUP in names(gene_set_groups)) {
        CATEGORY <- gene_set_groups[[GROUP]][1]
        PREFIX   <- gene_set_groups[[GROUP]][2]
        load(file.path(FGSEA_OUT, SAMPLE, paste0('set_', CATEGORY, '.Rdata')))
        sgnf <- fgseaRes %>% filter(padj < 0.05 & startsWith(pathway, PREFIX))
        # top_n keeps ties, as in the original notebook
        top <- c(top_n(sgnf, n = 10, wt = NES)$pathway, top_n(sgnf, n = -10, wt = NES)$pathway)
        if (length(top) > 0)
            candidates <- bind_rows(candidates, tibble(gene_set_group = GROUP, pathway = top, Sample = SAMPLE))
    }
}

candidates <- candidates %>%
    group_by(gene_set_group, pathway) %>%
    summarise(n_samples_top = n_distinct(Sample), .groups = 'drop') %>%
    mutate(in_fig3a = pathway %in% selected$pathway) %>%
    arrange(gene_set_group, desc(n_samples_top))

print(table(candidates$gene_set_group))
stopifnot(all(selected$pathway %in% candidates$pathway))
write.csv(candidates, OUT_CSV, row.names = FALSE)
