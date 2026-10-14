# Genes recurrently up-regulated in TFHL tumor cells vs normal TFH cells, across the 18 tumor samples.
# Concise version of ntbk_stat_deg.ipynb of the original analysis.
#
# A gene is significant in a sample when p_val_adj < 0.05, avg_log2FC > 1, pct.1 > 0.3 (expressed in
# > 30% of tumor cells) and pct.2 < 0.2 (expressed in < 20% of TFH cells).
#
# Usage: Rscript recurrent_deg.R <de_dir> [out_dir]
#   <de_dir>/<group>.csv from 02_differential_expression
# Output: <out_dir>/deg_sgnf_TFHL_vs_TFH.csv       significant genes, one row per gene x sample
#         <out_dir>/recurrent_deg_TFHL_vs_TFH.csv  genes significant in >= 2 samples, with the count

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
DE_DIR  <- args[1]
OUT_DIR <- if (length(args) >= 2) args[2] else '.'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df_deg <- NULL
for (f in sort(list.files(DE_DIR, pattern = '\\.csv$', full.names = TRUE))) {
    df <- read.csv(f, row.names = 1)
    df$gene <- rownames(df)
    df$sample <- sub('\\.csv$', '', basename(f))
    df_deg <- bind_rows(df_deg, df)
}
rownames(df_deg) <- NULL
print(table(df_deg$sample))

df_deg_sgnf <- df_deg %>% filter(p_val_adj < 0.05 & avg_log2FC > 1 & pct.1 > 0.3 & pct.2 < 0.2)
write.csv(df_deg_sgnf, file.path(OUT_DIR, 'deg_sgnf_TFHL_vs_TFH.csv'), row.names = FALSE)

recurrent <- df_deg_sgnf %>%
    count(gene, name = 'num_samples_sgnf_de') %>%
    filter(num_samples_sgnf_de >= 2) %>%
    arrange(desc(num_samples_sgnf_de), gene)
print(table(recurrent$num_samples_sgnf_de))
write.csv(recurrent, file.path(OUT_DIR, 'recurrent_deg_TFHL_vs_TFH.csv'), row.names = FALSE)

# The dot plot shows the genes significant in >= 8 samples, minus four removed by hand
# (see 04_dotplot/selected_genes.txt). Check that the list matches these results.
selected <- read.delim(file.path(HERE, '..', '04_dotplot', 'selected_genes.txt'), comment.char = '#')
top <- recurrent %>% filter(num_samples_sgnf_de >= 8)
cat('Genes significant in >= 8 samples and not in the dot plot:',
    setdiff(top$gene, selected$gene), '\n')
chk <- inner_join(selected, top, by = 'gene')
stopifnot(nrow(chk) == nrow(selected), all(chk$num_samples == chk$num_samples_sgnf_de))
