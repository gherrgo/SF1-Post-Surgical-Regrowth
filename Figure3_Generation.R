###############################################################################################
# --------------------------------------------------------------------------------------------#
# Figure 3a-c Generation - 
# This script provides classification and heatmap construction applicable across three publicly 
# available datasets:  Silva-Junior et al., 2023; Neou et al., 2020; Belakhoua et al., 2025.
# --------------------------------------------------------------------------------------------#
###############################################################################################
rm(list = ls()) %>% gc()

bioc_pkgs <- c( 'caret', 'e1071', 'IlluminaHumanMethylationEPICanno.ilm10b5.hg38', 'ComplexHeatmap')
cran_pkgs <- c( 'dplyr', 'magrittr', 'readr')

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

# ---- This syntax file will assume download & appropriate unpacking of publicly available sets
############################################################
# 1. Necessary materials acquisition
############################################################
load('./SF1_Classifier.Rdata') # Loads RF.obj into your environment - available through Mendeley DOI
dmp_sets <- read_excel('./Supplementary_Table_S2.xlsx')

############################################################
# 2. Publicly available data preprocessing & application of
# classifier. Processing standards chosen according to
# publication; end goal is a beta_matrix (e.g., silva_matrix).
############################################################

silva_prediction <- predict(RF.obj[[1]], t(beta_silva), type = 'prob', na.action = 'na.omit')
silva_prediction$Kcluster_5 <- colnames(silva_prediction)[max.col(silva_prediction, ties.method = "first")] # Assigning labels to the highest score per sample

############################################################
# 3. Heatmap construction - using previously generated DMPs
# to generate a heatmap across the newly assigned clusters
############################################################

# ----- Row Annotation construction
row.pt <- data.frame(dmp_sets$Memberships)
colnames(row.pt) <- c('Supervised Analysis')
colors <- list('Supervised Analysis' = c("k1" = "red", "k2" = "orange",
                                         "k3" = "yellow", 'k4' = 'green', 
                                         'k5' = 'blue'))
rowAnn.pt <- HeatmapAnnotation(df = row.pt, which = 'row', col = colors, 
                               gap = unit(0, 'mm'), show_annotation_name = FALSE, border = TRUE)

# ----- Column Annotation construction
col.pt <- data.frame(silva_prediction$Kcluster_5) # Here it is possible to add additional clinicopathological features
colnames(col.pt) <- c('Sample Assignment')
colors <- list('Sample Assignment' = c("k1" = "red", "k2" = "orange",
                                       "k3" = "yellow", 'k4' = 'green', 
                                       'k5' = 'blue'))
colAnn.pt <- HeatmapAnnotation(df = col.pt, which = 'col', col = colors,
                               gap = unit(0, 'mm'), show_annotation_name = FALSE, border = TRUE)

heatmap.pt = Heatmap(data.matrix(beta_silva),
                     name = "Methylation b-values",
                     col = matlab::jet.colors(200),
                     show_row_names = FALSE, cluster_rows = TRUE, 
                     cluster_columns = TRUE, column_split = col.pt$`Sample Assignment`,
                     show_row_dend = FALSE, show_column_names = FALSE,
                     top_annotation = colAnn.pt, show_heatmap_legend = FALSE,
                     right_annotation = rowAnn.pt, row_split = dmp_sets$Memberships,
                     row_title = "K-Cluster Specific DMPs - Silva-Junior et al.", 
                     row_title_side = "left", column_names_gp = gpar(fontsize = 10),
                     column_title_gp = gpar(fontsize = 10), show_column_dend = TRUE)
pdf('./Figures/Figure3a_SilvaHeatmap.pdf', height = 9, width = 11)
draw(heatmap.pt)
dev.off()