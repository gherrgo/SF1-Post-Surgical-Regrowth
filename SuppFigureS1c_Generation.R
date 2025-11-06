###############################################################################################
# --------------------------------------------------------------------------------------------#
# Supplementary Figure S1c generation - 
# This script allows for the classification of Fresh Frozen and FFPE PitNET tumor tissue samples
# and generation of a heatmap across the aforementioned k-cluster related DMP sets.
# --------------------------------------------------------------------------------------------#
###############################################################################################
rm(list = ls()) %>% gc()

bioc_pkgs <- c('caret', 'e1071', 'IlluminaHumanMethylationEPICanno.ilm10b5.hg38', 'minfi', 'magrittr',
               'ComplexHeatmap')
cran_pkgs <- c('tidyr', 'dplyr', 'readr', 'tibble', 'tidyverse')

f (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

############################################################
# 1. Necessary materials acquisition
############################################################

load('./Materials/Moller_KMeans_RandomForest_Classifier.RData') # Loads RF.obj into your environment - available through Mendeley DOI
dmp_sets <- read_excel('./Materials/Supplementary_File_S2.xlsx') #DMP information
load('./Materials/Unmasked_619k_Combat_Methyl.RData') # Loading the methylation matrix

pd <- read.xlsx('./FreshvsFFPE.xlsx', sheetName = 'Ark1')
pd <- pd[!is.na(pd$Sentrix_ID),]
pd$IDAT <- paste0(pd$Sentrix_ID, '_', pd$Sentrix_Position)
table(pd$IDAT %in% colnames(beta))
pd <- pd[pd$IDAT %in% colnames(beta),]
beta_freshvsffpe <- beta[, pd$IDAT]

rm(list = ls())
lapply(c('caret', 'e1071', 'minfi', 'IlluminaHumanMethylationEPICanno.ilm10b5.hg38', 
         'dplyr', 'magrittr', 'readr'), library, character.only = TRUE)

load('./Results/RF_Results_Fixed.RData') # Load your RF classifier
load('./Results/DMP_Information.RData') # Load your DMP sets
load('./Materials/Unmasked_MethylMatrix_GQSamples.RDS')
pd <- read.xlsx('./FreshvsFFPE.xlsx', sheetName = 'Ark1')
pd <- pd[!is.na(pd$Sentrix_ID),]
pd$IDAT <- paste0(pd$Sentrix_ID, '_', pd$Sentrix_Position)
table(pd$IDAT %in% colnames(beta))
pd <- pd[pd$IDAT %in% colnames(beta),]
beta_freshvsffpe <- beta[, pd$IDAT]

#######################################################################
# 2. Classification of fresh/frozen & FFPE samples according to cluster
#######################################################################
freshvsffpe_prediction <- predict(RF.obj[[1]], t(beta_freshvsffpe), type = 'prob', na.action = 'na.omit')
freshvsffpe_prediction$Kcluster_5 <- colnames(freshvsffpe_prediction)[max.col(freshvsffpe_prediction, ties.method = "first")]

freshvsffpe_clinical <- data.frame('Sample_ID' = c(colnames(beta_freshvsffpe)),
                                   'Kcluster_5' = freshvsffpe_prediction)
freshvsffpe_clinical$IDAT <- freshvsffpe_clinical$Sample_ID
freshvsffpe_clinical=merge(freshvsffpe_clinical,pd, by='IDAT')

##############################################################################
# 3. Generation of heatmap across aforementioned DMP sets (k-cluster specific)
##############################################################################

# ----- Generating a row annotation
row.pt <- data.frame(dmp_sets$Memberships)
colnames(row.pt) <- c('Supervised Analysis')
colors <- list('Supervised Analysis' = c("k1" = "red", "k2" = "orange",
                                         "k3" = "yellow", 'k4' = 'green', 'k5' = 'blue'))
rowAnn.pt <- HeatmapAnnotation(df = row.pt, which = 'row', col = colors, 
                               gap = unit(0, 'mm'), show_annotation_name = FALSE, border = TRUE)

# ----- Generating a column annotation ; K-means cluster assignment, Patient ID, Preparation and Lineage
names(freshvsffpe_clinical)[8] <- 'Kcluster_5'
col.pt <- data.frame(freshvsffpe_clinical$Kcluster_5, freshvsffpe_clinical$ID, freshvsffpe_clinical$Preparation, 
                     freshvsffpe_clinical$Lineage)
colnames(col.pt) <- c('Sample Assignment', 'Patient ID', 'Preparation', 'Lineage')
colors <- list('Sample Assignment' = c("1" = "red", "2" = "orange",
                                       "3" = "yellow", '4' = 'green', '5' = 'blue'),
               'Patient ID'= c("264" = "purple", "265" = "pink", "266" = "cyan","267" = "magenta",
                               "268" = "brown", "269" = "gray", "270" = "gold",       
                               "271" = "lightblue", "272" = "limegreen", "274" = "salmon"),
               'Preparation'= c("FFPE" ="red", "Frozen" = "dodgerblue"),
               'Lineage'= c("PIT1" = "yellow", "SF1" ="purple", "TPIT" = "blue"))
colAnn.pt <- HeatmapAnnotation(df = col.pt, which = 'col', col = colors,
                               gap = unit(0, 'mm'), show_annotation_name = TRUE, border = TRUE)

# ----- Heatmap construction
heatmap.pt = Heatmap(data.matrix(beta_freshvsffpe[dmp_sets$Probe_ID, ]),
                     name = "Methylation beta-values",
                     col = matlab::jet.colors(200),
                     show_row_names = FALSE,
                     cluster_rows = TRUE, 
                     cluster_columns = TRUE,
                     column_split = col.pt$`Sample Assignment`,
                     show_row_dend = FALSE,
                     show_column_names = FALSE,
                     top_annotation = colAnn.pt,
                     show_heatmap_legend = FALSE,
                     right_annotation = rowAnn.pt,
                     row_split = dmp_sets$Memberships,
                     row_title = "K-Cluster Specific DMPs - External Data", 
                     row_title_side = "left",
                     column_names_gp = gpar(fontsize = 10),
                     column_title_gp = gpar(fontsize = 10),
                     show_column_dend = TRUE)
draw(heatmap.pt)

pdf('./Figures/Supplemental_Figure_S1.pdf', height = 9, width = 11)
draw(heatmap.pt)
dev.off()
