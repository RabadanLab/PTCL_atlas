# Collect p-values and means of the selected ligand-receptor pairs between the
# T-cell population of interest and each TME cell type, from all CellPhoneDB runs.
# Concise version of ntbk_polished_dotplot_per_aitl.ipynb (data-assembly part).
#
# Usage: Rscript build_df_plot.R <cpdb_output_dir> [out.csv]
#   <cpdb_output_dir> must contain normal/, AITL/<sample>/ and NOS/<sample>/,
#   each holding pvalues.txt and means.txt (see 02_run_cellphonedb).

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(TRUE)
CPDB_OUT <- if (length(args) >= 1) args[1] else 'cpdb_output'
OUT_CSV  <- if (length(args) >= 2) args[2] else 'df_plot.csv'
HERE <- dirname(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE)))

tmes <- c('macrophage', 'cDC', 'endothelium', 'stromal', 'pDC', 'monocyte', 'B.cell', 'plasma',
          'CD4', 'CD8', 'Treg', 'GCB')

# CellPhoneDB stores each pair in one orientation only (partner_a_partner_b), so every
# selected pair has a "mirror" spelling; line i of both files refer to the same interaction.
interacting_pairs        <- read.table(file.path(HERE, 'selected_pairs.txt'))[, 1]
interacting_pairs_mirror <- read.table(file.path(HERE, 'selected_pairs_mirror.txt'))[, 1]

aitl_samples <- c('sn01', 'sn02', 'sn03', 'sn05', 'sn06', 'sn07', 'sn12', 'sn13', 'sn14',
                  'sn29', 'sn31', 'sn35', 'sn36', 'sn39', 'sn43', 'sn45', 'sn46', 'sn49')
nos_samples  <- c('sn30', 'sn32', 'sn33', 'sn34', 'sn37', 'sn40', 'sn41', 'sn44', 'sn48')

# Extract <focal cell type> <-> <TME> values for the selected pairs from one CellPhoneDB run.
#   focal:      cell-type label of the T-cell population in the CellPhoneDB meta file
#   label:      name of this column in the output (e.g. TFHL_sn01)
#   keep_missing: TRUE keeps NA rows when a pair/cell-type combination is absent
#                 (used for the normal reference, as in the original notebook)
extract_pairs <- function(run_dir, focal, label, keep_missing = FALSE) {
    df_p <- read.delim(file.path(run_dir, 'pvalues.txt'))
    df_m <- read.delim(file.path(run_dir, 'means.txt'))
    stopifnot(all(df_p$id_cp_interaction == df_m$id_cp_interaction))

    keep <- df_p$interacting_pair %in% c(interacting_pairs, interacting_pairs_mirror)
    df_p <- df_p[keep, ]; df_m <- df_m[keep, ]
    rownames(df_p) <- df_p$interacting_pair
    rownames(df_m) <- df_m$interacting_pair

    rows <- list()
    for (TME in tmes) {
        for (i in seq_along(interacting_pairs)) {
            pair <- interacting_pairs[i]
            col  <- paste0(focal, '.', TME)
            if (!pair %in% rownames(df_p)) {
                pair <- interacting_pairs_mirror[i]
                col  <- paste0(TME, '.', focal)
            }
            if (keep_missing || (pair %in% rownames(df_p) && col %in% colnames(df_p))) {
                rows[[length(rows) + 1]] <- data.frame(
                    interacting_pair = interacting_pairs[i], celltype = label, tme = TME,
                    p = df_p[pair, col], m = df_m[pair, col])
            }
        }
    }
    bind_rows(rows)
}

df_plot <- bind_rows(
    lapply(nos_samples, function(s)
        extract_pairs(file.path(CPDB_OUT, 'NOS', s), paste0('NOS_', s), paste0('NOS_', s))),
    extract_pairs(file.path(CPDB_OUT, 'normal'), 'CD4', 'CD4', keep_missing = TRUE),
    lapply(aitl_samples, function(s)
        extract_pairs(file.path(CPDB_OUT, 'AITL', s), 'T_tumor', paste0('TFHL_', s))),
    extract_pairs(file.path(CPDB_OUT, 'normal'), 'TFH', 'TFH', keep_missing = TRUE)
)

# p = 0 (no permutation exceeded) is floored at 1/iterations for -log10 display
df_plot$p[df_plot$p == 0] <- 0.001
df_plot$diagnosis <- sapply(strsplit(df_plot$celltype, '_'), `[`, 1)  # TFHL, NOS, CD4, TFH

write.csv(df_plot, OUT_CSV)
