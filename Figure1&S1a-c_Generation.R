###############################################################################################
# --------------------------------------------------------------------------------------------#
# Figure 1 & S1a-c Generation - 
# This script preprocesses methylation IDATs, performs batch correction,
# runs unsupervised k-means clustering, generates the ideal # of clusters (with user feedback)
# completes supervised analyses between clusters, and constructs both heatmaps and volcano plots
# --------------------------------------------------------------------------------------------#
###############################################################################################

rm(list = ls()) 
# ------------ Necessary package Installation & Loading
bioc_pkgs <- c(
  'minfi', 'IlluminaHumanMethylationEPICanno.ilm10b5.hg38', 'IlluminaHumanMethylationEPICmanifest',
  'ConsensusClusterPlus', 'sva', 'ComplexHeatmap', 'rtracklayer', 'annotatr'
)
cran_pkgs <- c(
  'dplyr', 'magrittr', 'readr', 'xlsx', 'reshape2', 'ggplot2', 'cowplot', 'stats'
)

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

############################################################
# 1. Preprocessing and Quality Control
############################################################

# 1. Preprocessing of samples - wd is set to folder with IDATs and clinical dataset
# GENERATING RGSET : MINFI
RGset <- read.metharray.exp(base = 'RawData/', force= T, recursive = T)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

############################################################
# 2. Probe Detection Filtering
############################################################

# DETECTION P-VALUE FILTERING OF CPG SITES 
detP <- detectionP(RGset, type = "m+u")
detP_failed  <- detP > 0.01
detP_remove <- names(which(rowMeans(detP_failed) > 0.5, TRUE))
print(length(detP_remove))

Mset <- preprocessIllumina(RGset)

# Probe removal as determined from p.detection
probes_removed <- match(detP_remove,rownames(Mset)) %>% unique %>% na.omit
Mset <- Mset[-probes_removed,]
dim(Mset) # This will print out the number of rows (sites) vs number of columns (samples)
qc <- getQC(Mset)

# Optional: Figure directory
dir.create('./Figures')
pdf('./Figures/Supplementary_Figure_1a.pdf', height = 6, width = 6)
plotQC(qc, badSampleCutoff = 10.5) # QC plotting
dev.off()

# FLAGGED SAMPLE : # 142 :: Visualize beta-distributions 
ggbg <- function() {
  points(0, 0, pch=16, cex=1e6, col="grey90")
  grid(col="white", lty=1)
}

par(mar=c(4,4,3,2), mgp=c(2.5,1,0), 
    cex.main=1.5, font.main="1", 
    fg="#6b6b6b", col.main="#4b4b4b")
pData(RGset)$Flagging <- c(rep('GQ', length.out = 141), 'BQ', rep('GQ', length.out = 50))
densityPlot(RGset, sampGroups = pData(RGset)$Flagging,
            main="Beta density plot", 
            xlab="Beta values", 
            panel.first=ggbg())

### Need to exclude this sample :: 142
Mset <- Mset[, -c(142)]

# SNP Probe removal :: (already completed on published data)
# drops the probes containing a SNP at the CpG interrogation and/or at the single nucleotide extension, for any minor allele frequency
 GRset <- mapToGenome(Mset)
# GRset <- addSnpInfo(GRset)
# GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)

# GENERATION OF BETA FROM RESULTING DATA - here we will use RGset, as SNPs are removed across published .idats 
beta <- data.frame(getBeta(GRset, type = 'EPIC'))
colnames(beta) <- sub('X', '', colnames(beta)) # Subbing the X out of the column names so we can map it exactly to the clinical data/manifest

# Reading in clinical data - we only want samples from the 'Main Cohort
pd <- read_excel("/media/data1/grayson/PeLB/Morten/Final_Syntax/Supplementary_File_S1.xlsx", 
                 sheet = "Main Cohort") # Supplementary File S1
beta <- beta[, pd$IDAT] # Subsetting the methylation beta based on samples present within the desired cohort

### ------- Unmasking procedure: quality control 
# NECESSARY OPERATORS:
'%notin%' <- Negate("%in%")
row_to_column <- function(df){
  var = colnames(df)[1]
  df <- df %>% remove_rownames() %>% column_to_rownames(var = var)
}
### Unmasking procedures and NA-value removal of entire genome _____
## Function built with consideration of EPIC hg38 Manifest
unmask <- function(df){
  #Loading EPIC manifest :: available for download at https://zwdzwd.github.io/InfiniumAnnotation 
  probe.anno <- read.table("EPIC.hg38.manifest.tsv.gz", sep = "\t", header = T)
  #Removal of masked probes
  probe.anno <- probe.anno[!probe.anno$MASK_general & probe.anno$CpG_chrm %in% paste0("chr", 1:22),]
  # Subsetting of probes found after unmasking
  df <- df[rownames(df) %in% probe.anno$probeID,]
  #Removal of any NA values
  df <- df[rowSums(is.na(df)) == 0,]
}

beta <- unmask(beta)

# Optional: Materials directory
dir.create('./Materials')
save(beta, file= './Materials/Unmasked_MethylMatrix_GQSamples.RData')

############################################################
# 3. Batch Correction (ComBat)
############################################################

# ---- Generation of batch clinical data
pd_batch <- pd[, c('IDAT', 'Batch')]
rownames(pd_batch) <- pd_batch$IDAT
colnames(pd_batch) <- c('sample_id', 'batch_id')

beta_combat <- ComBat(as.matrix(beta), batch = pd_batch$batch.id)
save(beta_combat, file = './Materials/Unmasked_619k_Combat_Methyl.RData')

# ----- Visualization of non-normalized and normalized data clustering through Principal Components
# ----- Non-normalized data: Supplementary Figure S1b
assay.pca <- prcomp(t(beta))
df_out <- as.data.frame(assay.pca$x)
df_out <- merge(df_out, pd[, c('IDAT', 'Batch', 'Lineage')], by.x = 0, by.y = 1)
df_out <- df_out %>% remove_rownames() %>% 
  column_to_rownames(var = 'Row.names')

# Custom ggplot2 plotting theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"),
               axis.ticks = element_line(color = "black"), plot.margin = unit(c(1,1,1,1), "line"))
# Percentage of variance delineation
percentage <- round(assay.pca$sdev / sum(assay.pca$sdev) *100, 2)
percentage <- paste(colnames(df_out3), "(", paste(as.character(percentage), "%", ")", sep = ""))

p1 <- ggplot(df_out3, aes(x=PC1, y=PC2, color = Batch)) + 
  geom_point(shape = 17, size = 4) + 
  theme_bw() + xlab(percentage[1]) + 
  ylab(percentage[2])+scale_color_manual(values = c("pink", "green"))

# ----- Normalized data: Supplementary Figure S1c
assay.pca <- prcomp(t(beta_combat))
df_normal <- as.data.frame(assay.pca$x)
df_normal <- merge(df_normal, pd[, c('IDAT', 'Batch', 'Lineage')], by.x = 0, by.y = 1)
df_normal <- df_normal %>% remove_rownames() %>% 
  column_to_rownames(var = 'Row.names')

# Percentage of variance delineation
percentage <- round(assay.pca$sdev / sum(assay.pca$sdev) *100, 2)
percentage <- paste(colnames(df_normal), "(", paste(as.character(percentage), "%", ")", sep = ""))

p2 <- ggplot(df_normal, aes(x=PC1, y=PC2, color = Batch)) + geom_point(shape = 17, size = 4) + 
  theme_bw() + xlab(percentage[1]) + ylab(percentage[2])+scale_color_manual(values = c("pink", "green"))

pdf('./Figures/SupplementaryFigure_S1cb.pdf', height = 5, width = 10)
p1+p2
dev.off()

# 3. Unsupervised grouping - using variant CpGs to define k-means driven clusters
beta.pt <- as.data.frame(t(beta))
variances <- apply(X=beta.pt, MARGIN = 2, FUN=var)
sorted <- sort(variances, decreasing = TRUE, index.return = TRUE)$ix[1:1000]
beta.pt <- beta[sorted, ]

############################################################
# 4. Unsupervised Consensus Clustering
############################################################

result_kmeans = ConsensusClusterPlus(as.matrix(beta.pt),maxK=12,reps=5000,pItem=0.7,pFeature=1, 
                                     title="./Figures/Unsupervised Clustering",clusterAlg="km",
                                     distance="euclidean",plot="pdf", 
                                     verbose = TRUE)

# ---- Optimal clustering # k = 5 - building k-cluster heatmap
consensus.5 <- result_kmeans[[5]]
col_fun = colorRamp2(c(0, .5, 1), c("white", "dodgerblue", "blue"))
ann.consensus <- data.frame(consensus.5$consensusClass)
names(ann.consensus)[1] <- "cluster"
consensus.colors <- list(cluster = c("1" = "red", "2" = "orange", "3" = "yellow", 
                                     "4" = "blue", '5' = 'green'))
colAnn.consensus <- HeatmapAnnotation(df = ann.consensus, which = 'col', col = consensus.colors, 
                                      gap = unit(0, 'mm'), show_annotation_name = TRUE, 
                                      gp = gpar(col="grey"), border = TRUE)

pdf('./Figures/K5_ClusteringHeatmap.pdf', height = 9, width = 11)
Heatmap(consensus.5$consensusMatrix, col = col_fun, column_split = consensus.5$consensusClass, 
        cluster_columns = TRUE, top_annotation = colAnn.consensus)
dev.off()

# Assigning results to the clinical data & saving
pd$Kcluster_5 <- consensus.5$consensusClass
write.xlsx(pd, file = './Materials/Clinical_KAssignment.RDS')

############################################################
# 5. Differentially Methylated Position (DMP) Analysis
############################################################
beta <- beta_combat

#######################################################
# ----------------------------------------------------#
#------------------------ k1 -------------------------#
# ----------------------------------------------------#
#######################################################

group1 <- c(pd$IDAT[pd$Kcluster_5 == 'k1']) 
group2 <- c(pd$IDAT[pd$Kcluster_5 != 'k1'])

# ------ Wilcoxon function definition
my.wilcox.test.p.value <- function(...){
  obj <- try(wilcox.test(...), silent = TRUE)
  if (is(obj, "try error"))
    return(NA)
  else
    return(obj$p.value)
}
# ---- p-value grabbing per CpG 
p <- apply(beta, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ---- Assigning Wilcoxon results
p.adj <- p.adjust(p, method = "fdr")
beta.tmp <- beta # Here we will assign the results to a temporary matrix
beta.tmp$p.value.raw <- p
beta.tmp$p.value.adj <- p.adj
beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2

#Setting thresholds of diff.means and pvalues
## Setting thresholds
p_val_threshold <- 0.05 # Traditional 95% confidence
pos_diffmean_threshold <- 0.35 
neg_diffmean_threshold <- -0.35

threshold_beta.tmp <- (beta.tmp$p.value.adj < p_val_threshold & (beta.tmp$diff.mean > pos_diffmean_threshold | beta.tmp$diff.mean < neg_diffmean_threshold))
length(which(threshold_beta.tmp)) # length of CpGs which hold significance
beta.tmp$threshold <- threshold_beta.tmp
beta.tmp$p.value.adj <- as.numeric(beta.tmp$p.value.adj)

# ---- Custom ggplot2 volcano plot theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), 
               axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
               plot.margin = unit(c(1,1,1,1), "line"))

############################################################
# 6. Volcano Plot Generation (k1)
############################################################

p1 <- ggplot(beta.tmp, aes(x = diff.mean, y = -1*log10(p.value.adj), color = threshold))+
  geom_point()+ theme + geom_vline(xintercept = c(neg_diffmean_threshold, pos_diffmean_threshold), linetype = "dotted")+
  geom_hline(yintercept = c(-1*log10(p_val_threshold)), linetype = "dotted")+
  ggtitle("KCluster 1 vs. Non-k1 Samples *Pituitary Adenoma Tissue*")+
  xlab(expression(paste("DNA Methylation difference (",beta.tmp,"-values)")))+
  ylab(expression(paste("-lo", g[10], "(adj p-values)")))+
  scale_fill_discrete(breaks = c("TRUE", "FALSE"))

# ------ Isolating and saving the results : k1
beta.tmp <- beta.tmp %>% filter(
  threshold == "TRUE"
)
k1_results <- beta_tmp[, c('p.value.adj', 'mean.group1', 'mean.group2', 'diff.mean')]
rownames(k1_results) <- rownames(beta_tmp)
dir.create("Results")
save(k1_results, file = './K1_Supervised_Statistics_DMPs.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k2 -------------------------#
# ----------------------------------------------------#
#######################################################

group1 <- c(pd$IDAT[pd$Kcluster_5 == 'k2']) 
group2 <- c(pd$IDAT[pd$Kcluster_5 != 'k2'])

# ---- p-value grabbing per CpG 
p <- apply(beta, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ---- Assigning Wilcoxon results
p.adj <- p.adjust(p, method = "fdr")
beta.tmp <- beta # Here we will assign the results to a temporary matrix
beta.tmp$p.value.raw <- p
beta.tmp$p.value.adj <- p.adj
beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2

## Setting thresholds
p_val_threshold <- 0.05 # Traditional 95% confidence
pos_diffmean_threshold <- 0.2
neg_diffmean_threshold <- -0.275

threshold_beta.tmp <- (beta.tmp$p.value.adj < p_val_threshold & (beta.tmp$diff.mean > pos_diffmean_threshold | beta.tmp$diff.mean < neg_diffmean_threshold))
length(which(threshold_beta.tmp)) # length of CpGs which hold significance
beta.tmp$threshold <- threshold_beta.tmp
beta.tmp$p.value.adj <- as.numeric(beta.tmp$p.value.adj)

# ---- Custom ggplot2 volcano plot theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), 
               axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
               plot.margin = unit(c(1,1,1,1), "line"))


############################################################
# 6. Volcano Plot Generation (k2)
############################################################

p2 <- ggplot(beta.tmp, aes(x = diff.mean, y = -1*log10(p.value.adj), color = threshold))+
  geom_point()+ theme + geom_vline(xintercept = c(neg_diffmean_threshold, pos_diffmean_threshold), linetype = "dotted")+
  geom_hline(yintercept = c(-1*log10(p_val_threshold)), linetype = "dotted")+
  ggtitle("KCluster 2 vs. Non-k2 Samples *Pituitary Adenoma Tissue*")+
  xlab(expression(paste("DNA Methylation difference (",beta.tmp,"-values)")))+
  ylab(expression(paste("-lo", g[10], "(adj p-values)")))+
  scale_fill_discrete(breaks = c("TRUE", "FALSE"))

# ------ Isolating and saving the results : k2
beta.tmp <- beta.tmp %>% filter(
  threshold == "TRUE"
)
k2_results <- beta_tmp[, c('p.value.adj', 'mean.group1', 'mean.group2', 'diff.mean')]
rownames(k2_results) <- rownames(beta_tmp)
save(k2_results, file = './k2_Supervised_Statistics_DMPs.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k3 -------------------------#
# ----------------------------------------------------#
#######################################################

group1 <- c(pd$IDAT[pd$Kcluster_5 == 'k3']) 
group2 <- c(pd$IDAT[pd$Kcluster_5 != 'k3'])

# ---- p-value grabbing per CpG 
p <- apply(beta, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ---- Assigning Wilcoxon results
p.adj <- p.adjust(p, method = "fdr")
beta.tmp <- beta # Here we will assign the results to a temporary matrix
beta.tmp$p.value.raw <- p
beta.tmp$p.value.adj <- p.adj
beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2

## Setting thresholds
p_val_threshold <- 0.05 # Traditional 95% confidence
pos_diffmean_threshold <- 0.275
neg_diffmean_threshold <- -0.20 

threshold_beta.tmp <- (beta.tmp$p.value.adj < p_val_threshold & (beta.tmp$diff.mean > pos_diffmean_threshold | beta.tmp$diff.mean < neg_diffmean_threshold))
length(which(threshold_beta.tmp)) # length of CpGs which hold significance
beta.tmp$threshold <- threshold_beta.tmp
beta.tmp$p.value.adj <- as.numeric(beta.tmp$p.value.adj)

# ---- Custom ggplot2 volcano plot theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), 
               axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
               plot.margin = unit(c(1,1,1,1), "line"))

############################################################
# 6. Volcano Plot Generation (k3)
############################################################

p3 <- ggplot(beta.tmp, aes(x = diff.mean, y = -1*log10(p.value.adj), color = threshold))+
  geom_point()+ theme + geom_vline(xintercept = c(neg_diffmean_threshold, pos_diffmean_threshold), linetype = "dotted")+
  geom_hline(yintercept = c(-1*log10(p_val_threshold)), linetype = "dotted")+
  ggtitle("KCluster 3 vs. Non-k3 Samples *Pituitary Adenoma Tissue*")+
  xlab(expression(paste("DNA Methylation difference (",beta.tmp,"-values)")))+
  ylab(expression(paste("-lo", g[10], "(adj p-values)")))+
  scale_fill_discrete(breaks = c("TRUE", "FALSE"))

# ------ Isolating and saving the results : k3
beta.tmp <- beta.tmp %>% filter(
  threshold == "TRUE"
)
k3_results <- beta_tmp[, c('p.value.adj', 'mean.group1', 'mean.group2', 'diff.mean')]
rownames(k3_results) <- rownames(beta_tmp)
save(k3_results, file = './k3_Supervised_Statistics_DMPs.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k4 -------------------------#
# ----------------------------------------------------#
#######################################################

group1 <- c(pd$IDAT[pd$Kcluster_5 == 'k4']) 
group2 <- c(pd$IDAT[pd$Kcluster_5 != 'k4'])

# ---- p-value grabbing per CpG 
p <- apply(beta, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ---- Assigning Wilcoxon results
p.adj <- p.adjust(p, method = "fdr")
beta.tmp <- beta # Here we will assign the results to a temporary matrix
beta.tmp$p.value.raw <- p
beta.tmp$p.value.adj <- p.adj
beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2

## Setting thresholds
p_val_threshold <- 0.05 # Traditional 95% confidence
pos_diffmean_threshold <- 0.5 
neg_diffmean_threshold <- -0.55 

threshold_beta.tmp <- (beta.tmp$p.value.adj < p_val_threshold & (beta.tmp$diff.mean > pos_diffmean_threshold | beta.tmp$diff.mean < neg_diffmean_threshold))
length(which(threshold_beta.tmp)) # length of CpGs which hold significance
beta.tmp$threshold <- threshold_beta.tmp
beta.tmp$p.value.adj <- as.numeric(beta.tmp$p.value.adj)

# ---- Custom ggplot2 volcano plot theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), 
               axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
               plot.margin = unit(c(1,1,1,1), "line"))


############################################################
# 6. Volcano Plot Generation (k4)
############################################################

p4 <- ggplot(beta.tmp, aes(x = diff.mean, y = -1*log10(p.value.adj), color = threshold))+
  geom_point()+ theme + geom_vline(xintercept = c(neg_diffmean_threshold, pos_diffmean_threshold), linetype = "dotted")+
  geom_hline(yintercept = c(-1*log10(p_val_threshold)), linetype = "dotted")+
  ggtitle("KCluster 4 vs. Non-k4 Samples *Pituitary Adenoma Tissue*")+
  xlab(expression(paste("DNA Methylation difference (",beta.tmp,"-values)")))+
  ylab(expression(paste("-lo", g[10], "(adj p-values)")))+
  scale_fill_discrete(breaks = c("TRUE", "FALSE"))

# ------ Isolating and saving the results : k4
beta.tmp <- beta.tmp %>% filter(
  threshold == "TRUE"
)
k4_results <- beta_tmp[, c('p.value.adj', 'mean.group1', 'mean.group2', 'diff.mean')]
rownames(k4_results) <- rownames(beta_tmp)
save(k4_results, file = './k4_Supervised_Statistics_DMPs.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k5 -------------------------#
# ----------------------------------------------------#
#######################################################

group1 <- c(pd$IDAT[pd$Kcluster_5 == 'k5']) 
group2 <- c(pd$IDAT[pd$Kcluster_5 != 'k5'])

# ---- p-value grabbing per CpG 
p <- apply(beta, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ---- Assigning Wilcoxon results
p.adj <- p.adjust(p, method = "fdr")
beta.tmp <- beta # Here we will assign the results to a temporary matrix
beta.tmp$p.value.raw <- p
beta.tmp$p.value.adj <- p.adj
beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2

## Setting thresholds
p_val_threshold <- 0.05 # Traditional 95% confidence
pos_diffmean_threshold <- 0.325 
neg_diffmean_threshold <- -0.275

threshold_beta.tmp <- (beta.tmp$p.value.adj < p_val_threshold & (beta.tmp$diff.mean > pos_diffmean_threshold | beta.tmp$diff.mean < neg_diffmean_threshold))
length(which(threshold_beta.tmp)) # length of CpGs which hold significance
beta.tmp$threshold <- threshold_beta.tmp
beta.tmp$p.value.adj <- as.numeric(beta.tmp$p.value.adj)

# ---- Custom ggplot2 volcano plot theme
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),
               strip.background = element_blank(), axis.text.x = element_text(color = "black"), 
               axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
               plot.margin = unit(c(1,1,1,1), "line"))

############################################################
# 6. Volcano Plot Generation
############################################################

p5 <- ggplot(beta.tmp, aes(x = diff.mean, y = -1*log10(p.value.adj), color = threshold))+
  geom_point()+ theme + geom_vline(xintercept = c(neg_diffmean_threshold, pos_diffmean_threshold), linetype = "dotted")+
  geom_hline(yintercept = c(-1*log10(p_val_threshold)), linetype = "dotted")+
  ggtitle("KCluster 5 vs. Non-k5 Samples *Pituitary Adenoma Tissue*")+
  xlab(expression(paste("DNA Methylation difference (",beta.tmp,"-values)")))+
  ylab(expression(paste("-lo", g[10], "(adj p-values)")))+
  scale_fill_discrete(breaks = c("TRUE", "FALSE"))

# ------ Isolating and saving the results : k5
beta.tmp <- beta.tmp %>% filter(
  threshold == "TRUE"
)
k5_results <- beta_tmp[, c('p.value.adj', 'mean.group1', 'mean.group2', 'diff.mean')]
rownames(k5_results) <- rownames(beta_tmp)
save(k5_results, file = './k5_Supervised_Statistics_DMPs.RData')

############################################################
# 7. Heatmap creation - DMP based
############################################################

# Subsetting those DMPs across the rows & creating a heatmap from them
beta.pt <- beta_combat[c(rownames(k1_results), rownames(k2_results), rownames(k3_results),
                         rownames(k4_results), rownames(k5_results)),]

# ------ Generation of column and row annotation 
pd_expanded <- data.frame(pd)
rownames(pd_expanded) <- pd_expanded$IDAT
pd_expanded <- pd_expanded[rownames(pd), ]

ann.pt <- data.frame(pd_expanded$Kcluster_5, pd_expanded$Kcluster_4, pd_expanded$InvasiveGrowth, pd_expanded$Primarysurgery, pd_expanded$Reinterventionagain)
names(ann.pt) <- c('KCluster5', 'KCluster4', 'Invasive Growth', 'Primary Surgery', 'Intervention Again')

colors <- list('KCluster5' = c("1" = "red", "2" = "orange", "3" = "yellow", "4" = "blue", '5' = 'green'),
               'KCluster4' = c("1" = "red", "2" = "orange", "3" = "yellow", "4" = "blue"),
               'Invasive Growth' = c('0' = 'white', '1' = 'black'),
               'Primary Surgery' = c('0' = 'white', '1' = 'black'),
               'Intervention Again' = c('0' = 'white', '1' = 'black'))
colAnn.pt <- HeatmapAnnotation(df = ann.pt, which = 'col', col = colors, 
                               gap = unit(0, 'mm'), show_annotation_name = TRUE, gp = gpar(col="grey"), border = TRUE)

dmp_sets <- data.frame('Probe_ID' = c(rownames(k1_results), rownames(k2_results), rownames(k3_results),
                                      rownames(k4_results), rownames(k5_results)),
                       'Memberships' = c(rep('k1', length(1:rownames(k1_results))),
                                         rep('k2', length(1:rownames(k2_results))),
                                         rep('k3', length(1:rownames(k3_results))),
                                         rep('k4', length(1:rownames(k4_results))),
                                         rep('k5', length(1:rownames(k5_results)))))

# ---- Row Annotation: CpG location information
gencode <- readGFF("./Materials/gencode.v38.annotation.gtf") # Comprehensive annotation download from https://www.gencodegenes.org/human/release_38.html 
gencode.s <- subset(gencode, type %in% "gene")  #select only "gene" category
gencode.g <- makeGRangesFromDataFrame(gencode.s, keep.extra.columns = T)

genes.promoters <- promoters(makeGRangesFromDataFrame(gencode.s, keep.extra.columns=T), upstream = 200, downstream = 200)
#get relevant probe annotations
probe.anno <- read.table("./Materials/EPIC.hg38.manifest.tsv.gz", sep = "\t", header = T) # Download from https://zwdzwd.github.io/InfiniumAnnotation 
probe.anno <- subset(probe.anno, Probe_ID %in% dmp_sets$Probe_ID)

#Creating a GRange object from probe annotations
annotations <- makeGRangesFromDataFrame(probe.anno, keep.extra.columns = TRUE, start.field = "CpG_beg", 
                                        end.field = "CpG_end", seqnames.field = "CpG_chrm", 
                                        strand.field = "probe_strand")

promoter_CpG <- findOverlaps(annotations, genes.promoters)
promoter_probes <- as.data.frame(subsetByOverlaps(annotations, genes.promoters, ignore.strand = F))
genebody_CpG <- findOverlaps(annotations, gencode.g)
genebody_probes <- as.data.frame(subsetByOverlaps(annotations,gencode.g, ignore.strand = F))

#CREATING ROW ANNOTATION FOR PROBES_______________________
dmp_sets$probefamily <- "Intergenic"
dmp_sets[dmp_sets$Probe_ID %in% genebody_probes$Probe_ID, "probefamily"] <- "Gene_body"
dmp_sets[dmp_sets$Probe_ID %in% promoter_probes$Probe_ID, "probefamily"] <- "Promoter"

row.pt <- data.frame(dmp_sets$probefamily)
colnames(row.pt) <- c('Probe family')

colors2 <- list('Probe family' = c("Gene_body" = "mistyrose", "Intergenic" = "limegreen",
                                   "Promoter" = "darkred"))
rowAnn.pt <- HeatmapAnnotation(df = row.pt, which = 'row', col = colors2, 
                               gap = unit(0, 'mm'), show_annotation_name = FALSE, border = TRUE)

### ENHANCER ANNOTATION
load("./Materials/genehancer_table_April2018.rda") # Custom compiled genehancer table, unavailable for download
#Creating a GRange object from probe annotations
annotations <- makeGRangesFromDataFrame(probe.anno, keep.extra.columns = TRUE, start.field = "CpG_beg", 
                                        end.field = "CpG_end", seqnames.field = "CpG_chrm",
                                        strand.field = "probe_strand")
#Creating a GRange object from Enhancer file
genehancer<- genehancer %>% filter(
  geneHancer.feature.name == "Enhancer"
)
enhancer <- makeGRangesFromDataFrame(genehancer, keep.extra.columns = TRUE, start.field = "geneHancer.start", end.field = "geneHancer.end", seqnames.field = "geneHancer.chrom",
                                     strand.field = "geneHancer.strand")
#Enhancer
enhancer_CpG <- findOverlaps(annotations, enhancer)
enhancer_probes <- as.data.frame(subsetByOverlaps(annotations, enhancer, ignore.strand = F))

dmp_sets$enhancer <- "No"
dmp_sets[dmp_sets$Probe_ID %in% enhancer_probes$Probe_ID, "enhancer"] <- "Yes"
row.pt2 <- data.frame(dmp_sets$enhancer)
colnames(row.pt2) <- c('Enhancer')

colors3 <- list(Enhancer = c("Yes" = "black", "No" = "white"))
rowAnn.pt2 <- HeatmapAnnotation(df = row.pt2, which = 'row', col = colors3, 
                                gap = unit(1, 'mm'), show_annotation_name = FALSE, border = TRUE)

#CpG LOCATION ANNOTATION ____________________________________________________________________________________________________________
annots = c('hg38_cpgs')
cpg_annotations = build_annotations(genome = 'hg38', annotations = annots)
dm_annotated = annotate_regions(
  regions = annotations,
  annotations = cpg_annotations, 
  ignore.strand = TRUE,
  quiet = FALSE
)
print(dm_annotated)
df_dm_annotated = data.frame(dm_annotated)

cpg_annotations = subset(df_dm_annotated, select = c(Probe_ID, annot.type))
cpg_annotations[cpg_annotations$annot.type == "hg38_cpg_inter", "CpG.Location"] <- "Open Seas"
cpg_annotations[cpg_annotations$annot.type == "hg38_cpg_shelves", "CpG.Location"] <- "Shelves"
cpg_annotations[cpg_annotations$annot.type == "hg38_cpg_shores", "CpG.Location"] <- "Shores"
cpg_annotations[cpg_annotations$annot.type == "hg38_cpg_islands", "CpG.Location"] <- "CpG Islands"
cpg_annotations = subset(cpg_annotations, select = c(Probe_ID, CpG.Location))
dmp_sets <- merge(dmp_sets, cpg_annotations, by = "Probe_ID", all.x = TRUE)
rownames(dmp_sets) <- dmp_sets$Probe_ID
dmp_sets <- dmp_sets[rownames(beta.pt),]

row.pt3 <- data.frame(dmp_sets$CpG.Location)
colnames(row.pt3) <- c('CpG Location')

colors4 <- list('CpG Location' = c("Open Seas" = "turquoise1", "Shelves" = "saddlebrown",
                                   "Shores" = "navajowhite", "CpG Islands" = "seagreen4"))
rowAnn.pt3 <- HeatmapAnnotation(df = row.pt3, which = 'row', col = colors4, 
                                gap = unit(1, 'mm'), show_annotation_name = FALSE, border = TRUE)

heatmap.pt = Heatmap(data.matrix(beta.pt),
                     name = "Methylation beta-values",
                     col = matlab::jet.colors(200),
                     show_row_names = FALSE,
                     cluster_rows = TRUE, 
                     cluster_columns = TRUE,
                     column_split = colAnn.pt$Kcluster_5,
                     show_row_dend = FALSE,
                     show_column_names = FALSE,
                     top_annotation = colAnn.pt,
                     show_heatmap_legend = FALSE,
                     right_annotation = rowAnn.pt,
                     row_split = dmp_sets$Memberships,
                     row_title = "K-Cluster Specific DMPs", row_title_side = "left",
                     column_names_gp = gpar(fontsize = 10),
                     column_title_gp = gpar(fontsize = 10),
                     show_column_dend = TRUE)
heatmap.pt = heatmap.pt + rowAnn.pt2 + rowAnn.pt3

pdf('./Figures/Figure1a_KMeans_DMP_Heatmap.pdf', height = 9, width = 11)
draw(heatmap.pt)
dev.off()

############################################################
# 8. Figure 1b compilation figure : volcano plots
############################################################

pdf('./Figures/Figure1b_DMPVolcanoPlots.pdf', height = 6, width = 12)
cowplot::plot_grid(p1, p2, p3, p4, p5, align = c('hv'), ncol = 5)
dev.off()