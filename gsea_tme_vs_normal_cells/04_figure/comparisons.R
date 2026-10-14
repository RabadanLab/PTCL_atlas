# The 8 comparisons (= DE groups), in the order used to collect results, and the two figure panels.
COMPARISONS <- c('Macrophage-AITL', 'Plasma-AITL', 'B-cell-AITL', 'cDC-AITL',
                 'Macrophage-NOS',  'Plasma-NOS',  'B-cell-NOS',  'cDC-NOS')
# Column order of each panel. The original notebook used ggplot's default alphabetical order
# (locale-dependent); it is fixed here.
PANEL_LEFT  <- c('B-cell-AITL', 'B-cell-NOS', 'Plasma-AITL', 'Plasma-NOS')
PANEL_RIGHT <- c('cDC-AITL', 'cDC-NOS', 'Macrophage-AITL', 'Macrophage-NOS')
