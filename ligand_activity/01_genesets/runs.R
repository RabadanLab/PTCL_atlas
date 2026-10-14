# The 6 NicheNet runs: 3 receiver cell types x 2 ecosystems (M1 = ecosystem 2, M2 = ecosystem 6).
# The gene set of interest of each run is the top 50 genes of one or two Spectra factors of the
# receiver cell type (columns of topgenes_per_factor.csv). For M1_to_stromal, the top 50 genes of
# two factors are combined.
# receiver: label in the `reanno` column of meta_annotation_Lin.csv
RUNS <- list(
    M1_to_macrophage  = list(receiver = 'macrophage',  factors = '35-X-Mac-X-35'),
    M2_to_macrophage  = list(receiver = 'macrophage',  factors = '37-X-Mac-X-37'),
    M1_to_endothelium = list(receiver = 'endothelium', factors = '73-X-endo-X-73'),
    M2_to_endothelium = list(receiver = 'endothelium', factors = '74-X-endo-X-74'),
    M1_to_stromal     = list(receiver = 'stromal',     factors = c('87-X-fibro-X-fibro_EGF_response',
                                                                   '93-X-fibro-X-93')),
    M2_to_stromal     = list(receiver = 'stromal',     factors = '91-X-fibro-X-fibro_TGFB1_response')
)
N_TOP_GENES <- 50
