###############################################################################################
# --------------------------------------------------------------------------------------------#
# Supplementary Figure S2 generation - 
# This script provides CNV oncoplot generation for the list of genes presented to be differentially
# altered among our SF-1 PitNET cohort.
# --------------------------------------------------------------------------------------------#
###############################################################################################
rm(list = ls()) %>% gc()

bioc_pkgs <- c('ExperimentHub', 'conumee2', 'plotly', 'minfi', 'ggplot2',
               'SummarizedExperiment', 'GenomicRanges', 'maftools', 'stringr', 'HGNChelper')
cran_pkgs <- c('tidyr', 'dplyr', 'readxl', 'tibble', 'tidyverse')

f (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

############################################################
# 1. Necessary materials acquisition
############################################################
# ----- Download EPIC Manifest: https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
# ----- Assuming materials are placed in prior 'Materials' folder
maniLocNew <- './Materials/EPIC_Manifest.csv'
Epic.mani.new <- read.csv(maniLocNew, sep = ",", header = T)

# ----- Download hg19 gencode (conumee compatible): https://www.gencodegenes.org/human/release_19.html
gencode.V19.Anno <- readGFF("./Materials/hg19.gff")

# ----- Necessary gene lists : Supplement File S2 - DMP derived genes, Supplementary File S3 - literature derived genes
geneList_DMPs <- read.delim('./Materials/Supplementary_File_S2.txt')
geneList_Lit <- read.delim('./Materials/Supplementary_File_S3.txt')

# ----- Clinicopathological information
pd <- read_excel('./Materials/Supplementary_File_S1.xlsx')
pit.idats <- pd$IDAT

# ----- Pituitary Control processing: taken from Herrgott, Grayson; Mosella, Maritza; Castro, Anavaleria (2025), “Raw DNA Methylation Data; Mosella et al., 2021 ”, Mendeley Data, V2, doi: 10.17632/5pzd2rg5ys.2 
# Assuming pituitary control IDATs have been downloaded (with Sample_Map.xlsx) and placed in 'Mosella_Data' folder
pitControls <- './Mosella_Data/'
mosella_samples <- read_excel('./Mosella_Data/Sample_Map.xlsx')
mosella_samples <- mosella_samples %>%
  filter(Sample.type == 'Nontumor') # 5 nontumor Pituitary controls

RGset.c <- read.metharray.exp(base = pitControls , recursive = F, force=T)
Mset.c <- preprocessIllumina(RGset.c)

# ----- Tumoral IDAT processing
RGset.q <- read.metharray.exp(base = 'RawData/', force= T, recursive = T)
Mset.q <- preprocessIllumina(RGset.q.sub)

############################################################
# 2. Conumee pipeline
############################################################
data.c <- CNV.load(Mset.c)
data.q <- CNV.load(Mset.q)

# ----- Setting genelist of interest: 'geneList_DMPs' or 'geneList_Lit'
geneList <- geneList.Lit

# ----- Generate gene list using aforementioned GENCODE
gene.rows <- length(geneList)
gene.columns <- 4
gene_df_Final <- as.data.frame(matrix(nrow = gene.rows, ncol = gene.columns))
rownames(gene_df_Final) = geneList
colnames(gene_df_Final) = c("Start","End", "Strand", "chr")

for (i in 1:length(rownames(gene_df_Final))) {
  currentGene <- as.character(rownames(gene_df_Final)[i])
  currentGene.test <- as.character(unique(gencode.V19.Anno$gene_name[gencode.V19.Anno$gene_name == currentGene]))
  
  if (identical(currentGene.test, currentGene) == 1) {
    
    gene_df_Final[i,1] <- min(unique(gencode.V19.Anno$start[gencode.V19.Anno$gene_name == currentGene]))  
    gene_df_Final[i,2] <- max(unique(gencode.V19.Anno$end[gencode.V19.Anno$gene_name == currentGene]))  
    gene_df_Final[i,3] <- as.character(unique(gencode.V19.Anno$strand[gencode.V19.Anno$gene_name == currentGene][1]))
    gene_df_Final[i,4] <- as.character(na.omit(unique(gencode.V19.Anno$seqid[gencode.V19.Anno$gene_name == currentGene]))[1])
    
    if (identical(currentGene.test, currentGene) == 1 & !is.na(gene_df_Final[i,4]) & !is.na(gene_df_Final[i,1]) & !is.na(gene_df_Final[i,2])) {
      
      grange.tmp <-GRanges(seqnames = as.character(gene_df_Final[i,4]), 
                           ranges = IRanges(start = as.numeric(gene_df_Final[i,1]), end = as.numeric(gene_df_Final[i,2])))
      
    } else {
      #skip 
    }
    
  } else {  
    
    next
  }
}

# ----- Drop if NA exists in either Start, End, or Chromosome
gene_df_Final.narm <- gene_df_Final %>% drop_na(c(Start,End,chr))
gene_df_Final.narm$name <- rownames(gene_df_Final.narm)

# ----- Removal of X/Y chromosome-linked genes
gene_df_Final.narm <- subset(gene_df_Final.narm, chr != "chrX")
gene_df_Final.narm <- subset(gene_df_Final.narm, chr != "chrY")

# ----- Removal of any remaining NA's
gene_df_Final.narm <- gene_df_Final.narm[!is.na(gene_df_Final.narm$arm), ]

# ----- Conversion to GRanges object
gene_df_Final.narm.gr <- makeGRangesFromDataFrame(gene_df_Final.narm, keep.extra.columns = T)

# ----- Defining regions to be excluded (e.g. polymorphic regions)
data(exclude_regions)

# ----- Creates CNV Annotation object
annoCNV <- CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06, 
                           array_type = "EPIC", exclude_regions = exclude_regions, # array_type = c("450k", "EPICv2") for analyzing EPICv2 (query) and 450k (controls) data
                           detail_regions = gene_df_Final.narm.gr, genome = "hg19", 
                           chrXY = FALSE)  

# ----- Normalization
genomSegCNV <- CNV.fit(data.q, data.c, annoCNV)
# ----- Genomic Binning
genomSegCNV <- CNV.bin(genomSegCNV)
# ----- Analysis of detail regions
genomSegCNV <- CNV.detail(genomSegCNV)
# ----- Segmentation
genomSegCNV <- CNV.segment(genomSegCNV)
# ----- Detection on focal CNVs
genomSegCNV <- CNV.focal(genomSegCNV)
names(genomSegCNV) <- sub('X', '', names(genomSegCNV)) # Subbing the X out of the column names so we can map it exactly to the clinical data/manifest

# ----- Write text outputs... depending on needs
segmentsCNV.segs <- CNV.write(genomSegCNV, what = "segments")
segmentsCNV.detail <- CNV.write(genomSegCNV, what = "detail")
segmentsCNV.bins <- CNV.write(genomSegCNV, what = "bins")
segmentsCNV.probes <- CNV.write(genomSegCNV, what = "probes")
segmentsCNV.gistic <- CNV.write(genomSegCNV, what = "gistic")

# ----- Storing CNV detailing object
segmentsCNV.detail <- CNV.write(genomSegCNV, what = "detail")

# ----- Final output : Literature derived CNVs
df.onco <- segmentsCNV.detail[complete.cases(segmentsCNV.detail), ] # Preserving only genes which have complete detailing across samples
lit_cnv <- ifelse(df.onco > 0.25, "Amp", 
                            ifelse(df.onco < -0.25, "Del", "Neutral"))
save(lit_cnv, file = './Materials/Literature_Derived_CNV.RData') # Optional saving
# ----- Process should be repeated for the DMP-derived genes as well. Named 'dmp_cnv'. Used in next step

####################################################################
# 3. Curation of relevant results into an oncoplot (Supp Figure S2)
####################################################################

# Necessary materials : dmp_cnv (cnv matrix for DMP related CNVs), lit_cnv (cnv matrix for Lit related CNVs), gene_anno (annotation for genes - Supplementary File S4)
total_cnv <- cbind(dmp_cnv, lit_cnv)
gene_anno <- read_excel('./Materials/Supplementary_File_S2.xlsx') # Resulting annotation curated internally

# ----- Gene Symbol updating: conversion to appropriate nomenclature (non-hg19) for chromosomal mapping w/ HGNChelper
genes <- colnames(total_cnv)
gene_fix <- checkGeneSymbols(genes)
gene_fix <- gene_fix[gene_fix$Approved == TRUE, ] # Only desire correctly labeled genes
total_cnv <- total_cnv[, gene_fix$x] # Subsetting of genes
colnames(total_cnv) <- gene_fix$Suggested.Symbol # Renaming according to appropriate symbols

# ----- Fixing nomenclature in gene_anno
gene_anno <- gene_anno[!duplicated(gene_anno$Gene), ] # Removing duplicated genes
colnames(gene_fix)[1] <- "Gene"
gene_anno <- gene_anno %>%
  left_join(gene_fix, by="Gene")
gene_anno$Gene <- gene_anno$Suggested.Symbol
gene_anno <- gene_anno[!is.na(gene_anno$Gene), ]
total_cnv <- total_cnv[, gene_anno$Gene] # Subsetting of total CNV matrix by the updated and available genes

# ---------------------------------------------------- #
# 3a. Conversion of cleaned CNV matrix into MAF format #
# ---------------------------------------------------- #
maf_like <- total_cnv %>%
  as.data.frame() %>%
  rownames_to_column("Tumor_Sample_Barcode") %>%
  pivot_longer(-Tumor_Sample_Barcode,
               names_to="Hugo_Symbol",
               values_to="Variant_Classification")

# --------------------------------------------------------- #
# 3b. Fetch coordinates from Ensembl for genes in CNV table #
# --------------------------------------------------------- #
mart <- useEnsembl("genes", dataset="hsapiens_gene_ensembl", mirror = 'www')

gene_coords <- getBM(
  filters="hgnc_symbol",
  attributes=c("hgnc_symbol", "chromosome_name",
               "start_position", "end_position", "band"),
  values=unique(maf_like$Hugo_Symbol),
  mart=mart
)
gene_coords <- gene_coords[gene_coords$band != '', ]
colnames(gene_coords) <- c("Hugo_Symbol","Chromosome",
                           "Start_Position","End_Position","band")
gene_coords$band <- substring(gene_coords$band, 1, 1)
gene_coords <- gene_coords %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Chromosome = ifelse(Chromosome %in% c(1:22,"X","Y"), Chromosome, "NA")
  )

# ----------------------------------------------------------------- #
# 3c. Merging into MAF-like format & preparation for MAF conversion #
# ----------------------------------------------------------------- #
maf_like2 <- maf_like %>%
  left_join(gene_coords, by="Hugo_Symbol") %>%
  drop_na(Chromosome, Start_Position, End_Position)

maf_fixed <- maf_like2 %>%
  mutate(
    Variant_Type = "CNV", Reference_Allele = "N", Tumor_Seq_Allele2 = "A"
  ) %>%
  dplyr::select(
    Hugo_Symbol, Chromosome, Start_Position, End_Position,
    Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode
  )

# ----- Creation of preliminary maf for CNV oncoplot
cnv_maf <- read.maf(
  maf = maf_fixed,
  verbose = TRUE, vc_nonSyn = c("Amp",
                                "Del")
)

# ------------------------------------------------------ #
# 3d. Generation of the oncoplot w/ necessary materials #
# ------------------------------------------------------ #

# ----- Setting the desired genes (and ensuring complete overlap)
gene_anno <- gene_anno %>% 
  dplyr::filter(Gene %in% cnv_maf@gene.summary$Hugo_Symbol)
desired_genes <- gene_anno$Gene

# ----- Adding desired clinical data to the MAF object
pd <- pd[, c('IDAT', 'Kcluster_5', 'Lineage')]
pd_ordered <- pd %>%
  filter(IDAT %in% cnv_maf@data$Tumor_Sample_Barcode) %>%
  arrange(factor(IDAT, levels = levels(cnv_maf@data$Tumor_Sample_Barcode)))
rownames(pd_ordered) <- pd_ordered$IDAT
colnames(pd_ordered)[1] <- 'Tumor_Sample_Barcode'
pd_ordered$Kcluster_5 <- as.character(pd_ordered$Kcluster_5)

# ----- Custom color scaling for oncoplot
cnv_colors <- c(
  "Amp" = "#E41A1C", "Del" = "#377EB8"
)

# ----- Final MAF creation
cnv_maf <- read.maf(
  maf = maf_fixed,
  verbose = TRUE, vc_nonSyn = c("Amp",
                                "Del"),
  clinicalData = pd_ordered 
)

# ----- Necessitating sample order by k-cluster
sample_order <- pd_ordered$Tumor_Sample_Barcode[order(pd_ordered$Kcluster_5, decreasing = FALSE)]

# ----- Oncoplot creation and saving
pdf('./Figures/Supplementary_Figure_S2.pdf', height = 11, width = 11)
oncoplot(
  maf = cnv_maf, genes = desired_genes, colors = cnv_colors,
  sampleOrder = sample_order, clinicalFeatures = c("Kcluster_5", "Lineage"),
  annotationColor = list(
    Kcluster_5 = c("1"="red","2"="orange","3"="yellow","4"="green","5"="blue"),
    Lineage = c('F' = '#000000', 'PIT1' = '#ffdd55', 'SF1' = '#a020f0', 'TPIT' = '#0000ff')),
  removeNonMutated = FALSE, drawRowBar = TRUE, drawColBar = TRUE,
  titleText = '84 differentially enriched & tumorigenesis involved genes across 117 PitNET samples'
)
dev.off()