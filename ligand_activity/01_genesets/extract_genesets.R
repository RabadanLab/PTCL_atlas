# Gene sets of interest of the 6 NicheNet runs, from the Spectra top-gene table.
# Replaces the `geneset_oi <- read.csv(...)` cells of the original
# {M1,M2}_to_{macrophage,endothelium,stromal}/ntbk_*_fixed_goi_rebak.ipynb notebooks.
#
# Usage: Rscript extract_genesets.R <topgenes_per_factor.csv> <out.csv>
#   topgenes_per_factor.csv: one column per Spectra factor, genes ranked by factor loading (200 rows)
# Output: CSV with columns run, receiver, factor, rank, gene. A gene shared by the two factors of
#         M1_to_stromal is listed once, under the first factor. The output of this script is
#         included as genesets_of_interest.csv.

HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))
source(file.path(HERE, 'runs.R'))

args <- commandArgs(TRUE)
# check.names = FALSE keeps the factor names as in the file (the notebooks used the
# read.csv-mangled names, e.g. X35.X.Mac.X.35)
topgenes <- read.csv(args[1], check.names = FALSE)

out <- do.call(rbind, lapply(names(RUNS), function(run) {
    r <- RUNS[[run]]
    df <- do.call(rbind, lapply(r$factors, function(f) {
        data.frame(run = run, receiver = r$receiver, factor = f,
                   rank = seq_len(N_TOP_GENES), gene = topgenes[[f]][seq_len(N_TOP_GENES)])
    }))
    df[!duplicated(df$gene), ]
}))
write.csv(out, args[2], row.names = FALSE, quote = FALSE)
