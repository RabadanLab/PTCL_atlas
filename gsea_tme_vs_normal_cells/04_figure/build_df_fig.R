# Collect the fgsea results of the figure's gene sets from the 8 TME-vs-normal comparisons.
# Concise version of the data-loading parts of ntbk_fgsea_tme.ipynb and ntbk_fgsea_tme_condensed.ipynb.
#
# Usage: Rscript build_df_fig.R <fgsea_output_dir> [out.csv]
#   <fgsea_output_dir>/<comparison>/set_{H,C2,C5}.Rdata from 03_fgsea
# Output: one row per gene set x comparison (fgsea columns except leadingEdge, plus Category, Sample).

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
FGSEA_OUT <- if (length(args) >= 1) args[1] else 'fgsea_output'
OUT_CSV   <- if (length(args) >= 2) args[2] else 'df_fgsea_tme_condensed.csv'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

source(file.path(HERE, 'comparisons.R'))
selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')

df <- NULL
for (COMPARISON in COMPARISONS) {
    for (CATEGORY in c('H', 'C2', 'C5')) {
        load(file.path(FGSEA_OUT, COMPARISON, paste0('set_', CATEGORY, '.Rdata')))
        df_tmp <- fgseaRes %>%
            filter(pathway %in% selected$pathway) %>%
            dplyr::select(-leadingEdge) %>%
            mutate(Category = CATEGORY, Sample = COMPARISON)
        df <- bind_rows(df, df_tmp)
    }
}

# A gene set can be absent from a comparison when fewer than 15 of its genes are ranked.
missing <- anti_join(expand_grid(pathway = selected$pathway, Sample = COMPARISONS), df,
                     by = c('pathway', 'Sample'))
if (nrow(missing) > 0) {
    message('Gene set / comparison combinations not tested by fgsea:')
    print(as.data.frame(missing))
}

write.csv(df, OUT_CSV, row.names = FALSE)
