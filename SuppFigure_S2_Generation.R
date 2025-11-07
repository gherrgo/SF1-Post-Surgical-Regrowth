###############################################################################################
# --------------------------------------------------------------------------------------------#
# Supplementary Figure S2a-e Generation - 
# This script provides differentially expressed gene (DEG) exploration as discovered through 
# application of the k-means classifier across Neou et al.'s matching methylome and transcriptome
# data collection.  Pipeline completes in identification of negatively correlated probe-gene pairs
# and generation of relevant scatterplots.
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
load('./Materials/Neou_Classification_Results.RData') # Results from the classification of Neou et al. through the methylation classifier

#####################################################################################
# 2. Publicly available data preprocessing  
# Assumes acquistion of raw transcription counts per sample (cts) and beta-matrix (e.g., beta_neou)
#####################################################################################
cts <- data.matrix(read.delim('./CountTable_RNAseq_PituitaryAdenoma.txt')) # Neou et al. samples count table
# -------- Exclusion of mitochondrial genome and ribosomal protein genes
exclusion_genes <- grep("^MT-|^RP", rownames(cts))
cts <- cts[-exclusion_genes, clinical_neou$`Source Name`]
identical(colnames(cts), clinical_neou$`Source Name`) # should be TRUE 

beta_neou <- beta_neou[, clinical_neou$`Source Name`]
identical(colnames(beta_neou), clinical_neou$`Source Name`) # should be TRUE 

#######################################################
# ----------------------------------------------------#
#------------------ DEG Derivation:  -----------------#
#---------------- k3 oriented analysis ---------------#
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

#####################################################################################
# 3. Linking of DEGs to regulatory CpGs - Enhancer and Promoter Regions (as defined)
#####################################################################################
probe.anno <- read.table("./Materials/EPIC.hg38.manifest.tsv.gz", sep = "\t", header = T) # Information on location can be found in Figure 1 File
probe.anno <- probe.anno[!probe.anno$MASK_general & probe.anno$CpG_chrm %in% paste0("chr", 1:22),]
probe.anno <- tidyr::separate_rows(probe.anno, gene, sep = ';', convert = FALSE)
probe.anno <- probe.anno[probe.anno$gene %in% rownames(k1_res05),]
potential_CpGs <- c(probe.anno$probeID)

# ------ Promoter mapping
gencode <- readGFF("./Materials/gencode.v38.annotation.gtf") # Information on location can be found in Figure 1 File
gencode.s <- subset(gencode, type %in% "gene")  #select only "gene" category
gencode.g <- makeGRangesFromDataFrame(gencode.s, keep.extra.columns = T)
# get promoters of these genes:
genes.promoters <- promoters(makeGRangesFromDataFrame(gencode.s, keep.extra.columns=T), upstream = 200, downstream = 200)

annotations <- makeGRangesFromDataFrame(probe.anno, keep.extra.columns = TRUE, start.field = "CpG_beg", end.field = "CpG_end", seqnames.field = "CpG_chrm",
                                        strand.field = "probe_strand")
promoter_CpG <- findOverlaps(annotations, genes.promoters)
promoter_probes <- as.data.frame(subsetByOverlaps(annotations, genes.promoters, ignore.strand = F))

# ------ Enhancer mapping
load("./Materials/genehancer_table_April2018.rda")  #Custom compiled genehancer table, unavailable for download
genehancer <- genehancer %>% filter(
  geneHancer.feature.name == "Enhancer"
)
enhancer <- makeGRangesFromDataFrame(genehancer, keep.extra.columns = TRUE, start.field = "geneHancer.start", end.field = "geneHancer.end", seqnames.field = "geneHancer.chrom",
                                     strand.field = "geneHancer.strand")
enhancer_CpG <- findOverlaps(annotations, enhancer)
enhancer_probes <- as.data.frame(subsetByOverlaps(annotations, enhancer, ignore.strand = F))

# ------ Constructing potential regulatory CpG list for exploration
potential_CpGs <- potential_CpGs[potential_CpGs %in% c(enhancer_probes$probeID, promoter_probes$probeID)]
potential_CpGs <- data.frame(CpG = c(potential_CpGs),
                             Type = c(ifelse(potential_CpGs %in% enhancer_probes$probeID, 'Enhancer', 'Promoter')))

#####################################################################################
# 4. Testing for differential methylation across Neou's methylation set
#####################################################################################
beta_neou <- beta_neou[potential_CpGs$CpG,]

group1 <- neou_clinical$`Source Name`[neou_clinical$KCluster_5 == '3'] # k3 samples for group1
group2 <- neou_clinical$`Source Name`[neou_clinical$KCluster_5 != '3'] # Non-k3 samples for group2

# ------ Our Non-parametric Wilcoxon function
my.wilcox.test.p.value <- function(...){
  obj <- try(wilcox.test(...), silent = TRUE)
  if (is(obj, "try error"))
    return(NA)
  else
    return(obj$p.value)
}
# ------ P-value derivation
p <- apply(beta_neou, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})

# ------ P-value adjustment and assigning wilcox results
p.adj <- p.adjust(p, method = "fdr") # p-value adjustment: False-discovery rate
potential_CpGs$p.value.raw <- p
potential_CpGs$p.value.adj <- p.adj
potential_CpGs$mean.group1 <- apply(beta_neou[, group1], 1, mean, na.rm = T)
potential_CpGs$mean.group2 <- apply(beta_neou[, group2], 1, mean, na.rm = T)
potential_CpGs$diff.mean <- potential_CpGs$mean.group1 - potential_CpGs$mean.group2

# ------ Thresholding of diff.means and pvalues : these are pliable
p_value_threshold <- 0.05
diff_mean_threshold <- 0.2

# ------ Optional filtering : here methylation beta-value difference is optional but ideal - you can remove if you find it further limiting your results
threshold_beta <- (poptential_CpGs$p.value.adj <= p_value_threshold & (abs(potential_CpGs$diff.mean) >= diff_mean_threshold))
length(which(threshold_beta))
potential_CpGs$threshold <- threshold_beta

# ------ Filtration of results according to thresholds
sig_results <- potential_CpGs[potential_CpGs$threshold == TRUE,]

#####################################################################################
# 5. Generation of subsequent scatter plots (DNA Methylation x Gene expression)
#####################################################################################
probe.anno2 <- probe.anno[probe.anno$probeID %in% sig_results$CpG, ]
probe.anno2$PGP <- paste0(probe.anno2$probeID, '_', probe.anno2$gene)

sig_results <- merge(probe.anno2[,c('probeID', 'gene', 'PGP')], sig_results, by = 1)
rownames(sig_results) <- sig_results$PGP

# ------ Scaling of count data for visualization
cts_neou <- data.frame(t(scale(t(as.matrix(cts))))) 
cts_neou <- cts_neou[sig_results$gene, c(group1, group2)] ## Necessary arranging of expression data
rownames(cts_neou) <- sig_results$PGP

beta_neou <- beta_neou[sig_results$CpG, c(group1, group2)] ## Necessary arranging of methylation data
rownames(beta_neou) <- sig_results$PGP

# ------ Convert beta_neou and cts_neou matrices into long format
beta_neou$PGP_ID <- rownames(beta_neou)
beta_long <- melt(beta_neou, varnames = c("PGP_ID", "Sample"), value.name = "beta_value")
cts_neou$PGP_ID <- rownames(cts_neou)
cts_long <- melt(cts_neou, varnames = c("PGP_ID", "Sample"), value.name = "expression_level")

# ------ Merge the two data frames based on 'PGP_ID' and 'Sample'
combined_df <- merge(beta_long, cts_long, by = c("PGP_ID", "variable"))

# ------ Subsetting pgps based off of negative correlations (Spearman testing)
cor_df <- combined_df %>%
  group_by(PGP_ID)%>% 
  summarize(corrrelation = cor(beta_value, expression_level, use = 'complete.obs', method = 'spearman'),
            p.value = cor.test(beta_value, expression_level, use = 'complete.obs', method = 'spearman')$p.value) %>%
  ungroup()
negative_cor_pgps <- cor_df$PGP_ID[cor_df$corrrelation < 0 & cor_df$p.value < 0.05]

# ------ Generating the scatterplot for the PGPs that we identify
combined_df  <- combined_df[combined_df$PGP_ID %in% negative_cor_pgps,]
combined_df$Group <- ifelse(combined_df$variable %in% group1, 'k3', 'Other')

pdf('./Figures/Supplementary_Figure_S2.pdf', height = 7, width = 10)
ggplot(combined_df, aes(x = beta_value, y = expression_level))+
  geom_point(aes(color = Group))+ 
  geom_smooth(method = 'lm', se = TRUE, color = 'orange') +
  labs(x = 'Methylation level (beta-value)', y = 'Expression Level (scaled)')+
  theme_bw() + stat_cor(method = 'spearman')+
  facet_wrap(~PGP_ID, scales = 'free') + scale_color_manual(values = c("orange", "black"))
dev.off()