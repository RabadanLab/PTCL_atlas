# Collect the fgsea results of the Fig. 3A gene sets from all 29 tumor-vs-CD4 comparisons.
# Concise version of the data-loading parts of ntbk_fgsea_tumor_vs_CD4.ipynb and ntbk_Fig3A.ipynb.
#
# Usage: Rscript build_df_fig3a.R <fgsea_output_dir> [out.csv]
#   <fgsea_output_dir>/<group>/set_{H,C2,C5}.Rdata from 03_fgsea
# Output: one row per gene set x sample (fgsea columns except leadingEdge, plus Category, Sample).

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
FGSEA_OUT <- if (length(args) >= 1) args[1] else 'fgsea_output'
OUT_CSV   <- if (length(args) >= 2) args[2] else 'df_Fig3A.csv'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

source(file.path(HERE, 'samples.R'))
selected <- read.delim(file.path(HERE, 'selected_pathways.txt'), comment.char = '#')

df <- NULL
for (SAMPLE in SAMPLES) {
    for (CATEGORY in c('H', 'C2', 'C5')) {
        load(file.path(FGSEA_OUT, SAMPLE, paste0('set_', CATEGORY, '.Rdata')))
        df_tmp <- fgseaRes %>%
            filter(pathway %in% selected$pathway) %>%
            dplyr::select(-leadingEdge) %>%
            mutate(Category = CATEGORY, Sample = SAMPLE)
        df <- bind_rows(df, df_tmp)
    }
}

# A gene set can be absent from a sample when fewer than 15 of its genes are ranked.
missing <- anti_join(expand_grid(pathway = selected$pathway, Sample = SAMPLES), df,
                     by = c('pathway', 'Sample'))
if (nrow(missing) > 0) {
    message('Gene set / sample combinations not tested by fgsea:')
    print(as.data.frame(missing))
}

write.csv(df, OUT_CSV, row.names = FALSE)
