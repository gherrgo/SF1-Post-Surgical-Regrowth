###############################################################################################
# --------------------------------------------------------------------------------------------#
# Figure 4 Generation - 
# This script provides differentially expressed gene (DEG) exploration as discovered through 
# application of the k-means classifier across Neou et al.'s matching methylome and transcriptome
# data collection.  Pipeline proceeds from SF-1 oriented DEG derivation to gene set enrichment
# (GSEA) completion.
# --------------------------------------------------------------------------------------------#
###############################################################################################
rm(list = ls()) %>% gc()

bioc_pkgs <- c('DESeq2', 'ggplot2', 'rtracklayer', 'clusterProfiler', 'pathview',
                'enrichplot', 'gprofiler2', 'org.Hs.eg.db', 'DOSE')
cran_pkgs <- c('tibble', 'dplyr')

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

# ------ This syntax file will assume download & appropriate unpacking of publicly available sets
############################################################
# 1. Necessary materials acquisition
############################################################
load('./SF1_Classifier.Rdata') # Loads RF.obj into your environment - available through Mendeley DOI

#####################################################################################
# 2. Publicly available data preprocessing & application of
# classifier. Processing standards chosen according to publication; end goal is
# a beta_matrix (e.g., beta_neou). Assumes acquistion of clinical data taken from the  
# publication (clinical_neou) and the raw transcription counts per sample (cts).
#####################################################################################
cts <- data.matrix(read.delim('./CountTable_RNAseq_PituitaryAdenoma.txt')) # Neou et al. samples count table
# -------- Exclusion of mitochondrial genome and ribosomal protein genes
exclusion_genes <- grep("^MT-|^RP", rownames(cts))
cts <- cts[-exclusion_genes,]

identical(colnames(beta_neou), clinical_neou$`Source Name`) # should be TRUE 
neou_prediction <- predict(RF.obj[[1]], t(beta_neou), type = 'prob', na.action = 'na.omit')
neou_clinical$Kcluster_5 <- colnames(neou_prediction)[max.col(neou_prediction, ties.method = "first")] # Assigning labels to the highest score per sample
save(neou_clinical, file = './Materials/Neou_Classification_Results.RData')

#####################################################################################
# 3. SF-1 only comparisons using DESeq2 across Neou et al raw expression data
# Includes 5 separate rounds of testing - one unique for each cluster
#####################################################################################
neou_clinical <- neou_clinical[neou_clinical$KCluster_5 != '4',] # Removal of 4th cluster (T-Pit & Pit-1)
cts <- cts[, neou_clinical$`Source Name`] # Ensuring only SF1 cluster presence
names(neou_clinical)[10] <- 'Lineage' # Some renaming

#######################################################
# ----------------------------------------------------#
#------------------------ k1 -------------------------#
# ----------------------------------------------------#
#######################################################
neou_clinical$k1_Others <- factor(ifelse(neou_clinical$KCluster == '1', 'k1', 'Others'), levels = c('Others', 'k1'))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = neou_clinical,
                              design = ~ k1_Others)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# ------ Filtration of genes by total attributed counts
smallestGroupSize <- min(table(neou_clinical$Lineage))
keep <- rowSums(counts(dds) >= 1) >= smallestGroupSize
table(keep) # # of genes being kept
dds <- dds[keep, ]

# ------ DESeq run::
dds <- DESeq(dds) # DESeq is run w/ default parameters
k1_res <- results(dds)
# ------ Fold-Change Shrinkage (optional):: Useful for visualization and ranking of genes in terms of fold-change between-groups
resLFC_k1 <- lfcShrink(dds, coef="k1_Others_k1_vs_Others", type="apeglm")

# ------ Refining results: Adjusted p-values & necessary filtration/saving
sum(k1_res$padj <= 0.05, na.rm=TRUE) #  Number of genes which are differentially expressed
k1_res05 <- results(dds, alpha=0.05)
summary(k1_res05)
k1_res05 <- k1_res[!is.na(k1_res$padj) & k1_res$padj <= 0.05,]
save(k1_res05, file = './Results/DEGs_k1_ONLYSF1.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k2 -------------------------#
# ----------------------------------------------------#
#######################################################
neou_clinical$k2_Others <- factor(ifelse(neou_clinical$KCluster == '2', 'k2', 'Others'), levels = c('Others', 'k2'))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = neou_clinical,
                              design = ~ k2_Others)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# ------ Filtration of genes by total attributed counts
smallestGroupSize <- min(table(neou_clinical$Lineage))
keep <- rowSums(counts(dds) >= 1) >= smallestGroupSize
table(keep) # # of genes being kept
dds <- dds[keep, ]

# ------ DESeq run::
dds <- DESeq(dds) # DESeq is run w/ default parameters
k2_res <- results(dds)
# ------ Fold-Change Shrinkage (optional):: Useful for visualization and ranking of genes in terms of fold-change between-groups
resLFC_k2 <- lfcShrink(dds, coef="k2_Others_k2_vs_Others", type="apeglm")

# ------ Refining results: Adjusted p-values & necessary filtration/saving
sum(k2_res$padj <= 0.05, na.rm=TRUE) #  Number of genes which are differentially expressed
k2_res05 <- results(dds, alpha=0.05)
summary(k2_res05)
k2_res05 <- k2_res[!is.na(k2_res$padj) & k2_res$padj <= 0.05,]
save(k2_res05, file = './Results/DEGs_k2_ONLYSF1.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k3 -------------------------#
# ----------------------------------------------------#
#######################################################
neou_clinical$k3_Others <- factor(ifelse(neou_clinical$KCluster_5 == '3', 'k3', 'Others'), levels = c('Others', 'k3'))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = neou_clinical,
                              design = ~ k3_Others)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# ------ Filtration of genes by total attributed counts
smallestGroupSize <- min(table(neou_clinical$Lineage))
keep <- rowSums(counts(dds) >= 1) >= smallestGroupSize
table(keep) # # of genes being kept
dds <- dds[keep, ]

# ------ DESeq run::
dds <- DESeq(dds) # DESeq is run w/ default parameters
k3_res <- results(dds)
# ------ Fold-Change Shrinkage (optional):: Useful for visualization and ranking of genes in terms of fold-change between-groups
resLFC_k3 <- lfcShrink(dds, coef="k3_Others_k3_vs_Others", type="apeglm")

# ------ Refining results: Adjusted p-values & necessary filtration/saving
sum(k3_res$padj <= 0.05, na.rm=TRUE) #  Number of genes which are differentially expressed
k3_res05 <- results(dds, alpha=0.05)
summary(k3_res05)
k3_res05 <- k3_res[!is.na(k3_res$padj) & k3_res$padj <= 0.05,]
save(k3_res05, file = './Results/DEGs_k3_ONLYSF1.RData')

#######################################################
# ----------------------------------------------------#
#------------------------ k5 -------------------------#
# ----------------------------------------------------#
#######################################################
neou_clinical$k5_Others <- factor(ifelse(neou_clinical$KCluster == '5', 'k5', 'Others'), levels = c('Others', 'k5'))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = neou_clinical,
                              design = ~ k5_Others)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# ------ Filtration of genes by total attributed counts
smallestGroupSize <- min(table(neou_clinical$Lineage))
keep <- rowSums(counts(dds) >= 1) >= smallestGroupSize
table(keep) # # of genes being kept
dds <- dds[keep, ]

# ------ DESeq run::
dds <- DESeq(dds) # DESeq is run w/ default parameters
k5_res <- results(dds)
# ------ Fold-Change Shrinkage (optional):: Useful for visualization and ranking of genes in terms of fold-change between-groups
resLFC_k5 <- lfcShrink(dds, coef="k5_Others_k5_vs_Others", type="apeglm")

# ------ Refining results: Adjusted p-values & necessary filtration/saving
sum(k5_res$padj <= 0.05, na.rm=TRUE) #  Number of genes which are differentially expressed
k5_res05 <- results(dds, alpha=0.05)
summary(k5_res05)
k5_res05 <- k5_res[!is.na(k5_res$padj) & k5_res$padj <= 0.05,]
save(k5_res05, file = './Results/DEGs_k5_ONLYSF1.RData')

#####################################################################################
# 4. GSEA on k3 relevant DEGs (derived from prior analysis) - cluster of interest.
#####################################################################################

# ------ Gene ranking according to fold-change
k3_DEGs <- data.frame(k3_res05)
gene_list <- c(k3_DEGs$log2FoldChange)
names(gene_list) <- rownames(k3_DEGs)
gene_list <- sort(gene_list, decreasing = TRUE)

# ------ Figure 4a : clusterProfiler barplot
names <- names(gene_list)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

gse <- gseGO(geneList=gene_list,  ont ="ALL", 
             keyType = "SYMBOL",  nPerm = 10000, 
             minGSSize = 3,  maxGSSize = 500, 
             pvalueCutoff = 0.05,  verbose = TRUE, 
             OrgDb = organism,  pAdjustMethod = "none") # No adjustment due to small sample size

edox <- setReadable(gse, 'org.Hs.eg.db', 'SYMBOL')
pdf('./Figures/Figure4a.pdf', height = 6, width = 9)
cnetplot(gse, categorySize = 'pvalue', foldChange = gene_list, circular = TRUE, colorEdge = TRUE)
dev.off()

# ------ Figure 4b : clusterProfiler barplot
names <- bitr(names,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = organism)
egobp <- enrichGO(
  gene = names[[2]], OrgDb    = 'org.Hs.eg.db',
  ont = "ALL", pAdjustMethod = "none",
  pvalueCutoff = 0.05, readable = TRUE)

pdf('./Figures/Figure4b.pdf', height = 6, width = 8)
ggplot(egobp[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))+theme_bw()
dev.off()
