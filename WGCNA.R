#=============================================================================#
# WGCNA.R                                                                     #
#	                                                                            #
# Start date: 11/01/2026      									                              #
#=============================================================================#

# Load packages ------------------------------------------------
library(rstudioapi)
library(tidyr)
library(dplyr)
library(missMethods)
library(ggplot2)

library(WGCNA)
allowWGCNAThreads()
library(magick)

library(RCy3)
library(readr)
library(readxl)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#-----------------------------------------------------------------------------#
# Data import (same as in DataPreparationExploratoryAnalysis.R)
#-----------------------------------------------------------------------------#
# Set working directory according to the location of this script
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load gene expression data and participant information (metadata)
# NOTE - *SVA-corrected* gene expression dataset (log2-transformed CPM)
gxData_all <- read.csv("Data/CPMS_SVA_corrected.csv", as.is = T, row.names = 1)
gxData_all <- log2(gxData_all + 1)  #Log2 transform
gxData_all <- impute_median(gxData_all, type = "columnwise")
metadata_all <- read.csv("Data/MAGNET_GX_2025/MAGNET_SampleData_18112022.csv", as.is = T, row.names = 1)

# Filter out only female participants who either have DCM or are healthy
metadata <- metadata_all |> 
  filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))

# Filter gxData_all to include only participants also in metadata
gxData <- gxData_all[, colnames(gxData_all) %in% rownames(metadata)]

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

#-----------------------------------------------------------------------------#
# Weighed Gene-Coexpression Network Analysis (WGCNA)
#-----------------------------------------------------------------------------#
# Basic quality check
gsg <- goodSamplesGenes(gxData, verbose = 3)

if (gsg$allOK) {
  cat("✓ All genes and samples passed quality checks!\n")
} else {
  cat("⚠ Some genes or samples failed quality checks.\n")
  # Print details if failed
  if (sum(!gsg$goodGenes) > 0) {
    cat("  - Genes failing:", sum(!gsg$goodGenes), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("  - Samples failing:", sum(!gsg$goodSamples), "\n")
  }
}

# Message should say: "✓ All genes and samples passed quality checks!"

# Transpose: WGCNA requires samples as rows, genes as columns
gxData_t <- as.data.frame(t(gxData))

# Check if samples are in the same order in metadata and gxData_t
all(rownames(metadata) == rownames(gxData_t)) 


# Encoding categorical variables ---------------------------------------------#
# Remove unnecessary columns
metadata$gender <- NULL
metadata$disease_race <- NULL
metadata$Library.Pool <- NULL # NOTE - Library.pool was removed!

# tissue source
metadata$tissue_source <- as.character(metadata$tissue_source)
metadata$tissue_source[metadata$tissue_source == "Cardiectomy"] <- 1

# etiology
metadata$etiology <- as.character(metadata$etiology)
metadata$etiology[metadata$etiology == "NF"] <- 0
metadata$etiology[metadata$etiology == "DCM"] <- 1

# race
metadata$race <- as.character(metadata$race)
metadata$race[metadata$race == "AA"] <- 0
metadata$race[metadata$race == "Caucasian"] <- 1

# afib
metadata$afib <- as.character(metadata$afib)
metadata$afib[metadata$afib == "No"] <- 0
metadata$afib[metadata$afib == "Yes"] <- 1

# VTVF
metadata$VTVF <- as.character(metadata$VTVF)
metadata$VTVF[metadata$VTVF == "No"] <- 0
metadata$VTVF[metadata$VTVF == "Yes"] <- 1

# Diabetes
metadata$Diabetes <- as.character(metadata$Diabetes)
metadata$Diabetes[metadata$Diabetes == "No"] <- 0
metadata$Diabetes[metadata$Diabetes == "Yes"] <- 1

# Hypertension
metadata$Hypertension <- as.character(metadata$Hypertension)
metadata$Hypertension[metadata$Hypertension == "No"] <- 0
metadata$Hypertension[metadata$Hypertension == "Yes"] <- 1

# Convert all columns to numeric (except tissue_source, race, and Library.pool)
lapply(list(metadata$etiology,
            metadata$age, 
            metadata$weight, 
            metadata$height, 
            metadata$hw, 
            metadata$lv_mass,
            metadata$afib,
            metadata$VTVF,
            metadata$Diabetes,
            metadata$Hypertension,
            metadata$LVEF,
            metadata$RIN,
            metadata$minexpr, 
            metadata$TIN.median.), as.numeric)

# Make all variables numeric
metadata$tissue_source <- NULL # gives problems when using mutate_all
metadata <- mutate_all(metadata, function(x) as.numeric(as.character(x)))

# Display trait summary
summary(metadata)

# Clean up memory
collectGarbage()


# Sample clustering and trait visualisation ----------------------------------#
# Hierarchical clustering of samples based on gene expression
sampleTree <- hclust(dist(gxData_t), method = "average")

meta.traits <- `colnames<-`(
  metadata[, c("etiology", "age", "race", "LVEF", "VTVF", "Hypertension", "Diabetes")],
  c("DCM", "Age", "Race", "LVEF", "VTVF", "Hypertension", "Diabetes"))
traitColors <- numbers2colors(meta.traits, signed = FALSE) # Convert those traits to colors

# Plot and save dendrogram with trait heatmap
dendrogram_path <- file.path("/Users/mikiverme/Desktop/DCM_diagnostics", "dendrogram_heatmap.png")
png(filename = dendrogram_path, width = 2400, height = 1800, res = 300)
plotDendroAndColors(
  sampleTree, 
  traitColors,
  groupLabels = names(meta.traits), 
  cex.dendroLabels = 0.5, 
  main = "Sample Dendrogram and Trait Heatmap"
)
dev.off()


# Network construction - Soft-thresholding -----------------------------------#
# Test different powers
powers <- seq(1, 15, by = 1)
sft <- pickSoftThreshold(gxData_t, powerVector = powers, verbose = 5) # Analyze network topology for each power

save(sft, file = "WGCNA-sft.RData")
sft$fitIndices # display result in table

# Test sft power mean connectivity
targetPower.test <- 5
meanConnectivity <- sft$fitIndices[               # Test of mean connectivity for power '4'
  sft$fitIndices[, "Power"] == targetPower.test,
  "mean.k."
]
meanConnectivity


# Visualise power selection and degree distribution --------

png(filename = "powerplots.png", width = 2400, height = 800, res = 300)

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(5, 4, 4, 2) + 0.1)

# Plot 1: Scale-free topology fit
plot(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (signed R²)",
  type = "n",
  main = "Scale Independence"
)
title(main = "A)", adj = 0, col.main = "grey50", font.main = 2, cex.main = 1.5)
text(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
  labels = powers,
  cex = 0.9,
  col = "red"
)
abline(h = 0.85, col = "red", lty = 2)

# Plot 2: Mean connectivity
plot(
  sft$fitIndices[, 1], 
  sft$fitIndices[, 5], 
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity", 
  type = "n", 
  main = "Mean Connectivity"
)
title(main = "B)", adj = 0, col.main = "grey50", font.main = 2, cex.main = 1.5)
text(
  sft$fitIndices[, 1], 
  sft$fitIndices[, 5], 
  labels = powers, 
  cex = 0.9,
  col = "red"
)

# Plot 3: Degree Distribution
hist(
  network.metrics$Degree,
  breaks = 30,
  main = "Degree Distribution",
  xlab = "Degree",
  ylab = "Number of Nodes"
)
title(main = "C)", adj = 0, col.main = "grey50", font.main = 2, cex.main = 1.5)

par(op)
dev.off()


# Build co expression network ------------------------------------------------#
rm(cor)

# Construct the network and identify modules
net <- blockwiseModules(
  gxData_t,
  power = 5,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "expTOM",
  verbose = 3,
  corType = "pearson",
  corOptions = list(use = "pairwise.complete.obs")
)

cat("  - Number of modules found:", length(unique(net$colors)), "\n")
# Number of modules found the first time (11/01/2026) --> 39

# Save network object
save(net, file = "WGCNA-net.RData")


# Module identification and visualisation ------------------------------------#

# Convert numeric module labels to colors
mergedColors <- labels2colors(net$colors)

# GENE BLOCK 1
# Plot gene dendrogram with module colors below per block
png("block1_dendrogram.png", width = 1200, height = 600, res = 150)                    
plotDendroAndColors(
  net$dendrograms[[1]], 
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors (block 1)"
)
dev.off()

# GENE BLOCK 2
png("block2_dendrogram.png", width = 1200, height = 600, res = 150)                    
plotDendroAndColors(
  net$dendrograms[[2]], 
  mergedColors[net$blockGenes[[2]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors (block 2)"
)
dev.off()

# GENE BLOCK 3
png("block3_dendrogram.png", width = 1200, height = 600, res = 150)                    
plotDendroAndColors(
  net$dendrograms[[3]], 
  mergedColors[net$blockGenes[[3]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors (block 3)"
)
dev.off()

# GENE BLOCKS TOGETHER (Dendrograms)
save_block <- function(block_num, label) {
  png(paste0("block", block_num, ".png"), width = 1200, height = 600, res = 150)
  
  plotDendroAndColors(
    net$dendrograms[[block_num]],
    mergedColors[net$blockGenes[[block_num]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = ""
  )
  
  # Draw label on top-left, in front of plot
  par(xpd = NA)
  text(0, 6, label, adj = c(0,1), col = "grey50", font = 2, cex = 1.5)
  
  dev.off()
}

save_block(1, "A)")
save_block(2, "B)")
save_block(3, "C)")


# Read images
b1 <- image_read("block1.png")
b2 <- image_read("block2.png")
b3 <- image_read("block3.png")

combined <- c(b1, b2, b3) %>% image_append(stack = TRUE)
image_write(combined, "combined_blocks.png")

# Store module information
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree1 <- net$dendrograms[[1]]
geneTree2 <- net$dendrograms[[2]]
geneTree3 <- net$dendrograms[[3]]

# Count genes per module
module.table <- table(moduleColors)
module.df <- as.data.frame(module.table)
colnames(module.df) <- c("Module", "Number of Genes")
module.df <- module.df[order(module.df$`Number of Genes`, decreasing = TRUE), ]

print(module.df)


# Module Trait Relationships ----------------------CORRECTED-------------------#

# Calculate module eigengenes (MEs) - first principal component of each module
nGenes <- ncol(gxData_t)
nSamples <- nrow(gxData_t)

MEs0 <- moduleEigengenes(gxData_t, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Select traits of interest explicitly
traits_of_interest <- metadata[, c("etiology", "age", "race", "LVEF", "VTVF", "Diabetes", "Hypertension")]  # adjust as needed

# Convert all traits to numeric
traits_of_interest <- mutate_all(traits_of_interest, function(x) as.numeric(as.character(x)))

# Correlate MEs with multiple traits
moduleTraitCor <- cor(MEs, traits_of_interest, use = "pairwise.complete.obs")
rownames(moduleTraitCor) <- colnames(MEs)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Extract correlations and p-values for etiology trait
dcMcor <- moduleTraitCor[, "etiology"]
dcMp <- moduleTraitPvalue[, "etiology"]

# Rank modules by strength of correlation with etiology
ranked.modules <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = dcMcor,
  Pvalue = dcMp
)
ranked.modules <- ranked.modules[order(-abs(ranked.modules$Correlation)), ]

print(ranked.modules)

# Subset the top 5 ranked modules (for later)
top5.modules <- ranked.modules$Module[!ranked.modules$Module %in% c("MEdarkgrey", "MEgrey")][1:5]


# Heatmap ---------------------------------------------------------------------#

# Create text matrix with correlations and p-values
moduleTraitCor_mat <- matrix(moduleTraitCor, ncol = 1)
colnames(moduleTraitCor_mat) <- "etiology"
rownames(moduleTraitCor_mat) <- names(moduleTraitCor)

moduleTraitPvalue_mat <- matrix(moduleTraitPvalue, ncol = 1)
colnames(moduleTraitPvalue_mat) <- "etiology"
rownames(moduleTraitPvalue_mat) <- names(moduleTraitPvalue)

moduleTraitCor_sub <- moduleTraitCor[top5.modules, , drop = FALSE]
moduleTraitPvalue_sub <- moduleTraitPvalue[top5.modules, , drop = FALSE]

textMatrix_sub <- matrix(
  paste(
    signif(moduleTraitCor_sub, 2),
    "\n(",
    signif(moduleTraitPvalue_sub, 1),
    ")",
    sep = ""
  ),
  nrow = nrow(moduleTraitCor_sub),
  ncol = ncol(moduleTraitCor_sub)
)

# Create heatmap
ySymbols_sub <- c("Blue", "Turquoise", "Brown", "Yellow", "Green")  # example manual names
yLabels_sub <- top5.modules
xLabels_sub <- c("DCM", "Age", "Race", "LVEF", "VTVF", "Hypertension", "Diabetes")  # match order of columns in moduleTraitCor_sub

png("top5_modules_heatmap.png", width = 2000, height = 1600, res = 300)
par(mar = c(8, 8.5, 3, 3))  # Adjust margins for your plot

labeledHeatmap(
  Matrix = moduleTraitCor_sub,
  xLabels = xLabels_sub,
  yLabels = yLabels_sub,
  ySymbols = ySymbols_sub,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_sub,
  setStdMargins = FALSE,
  cex.text = 1,
  cex.lab.y = 0.6,
  zlim = c(-1, 1),
  main = "Top 5 Module-Trait Relationships"
)

dev.off()


# Gene module membership and significance ------------------------------------#

# Define top modules (from  ranked.modules output)
modules <- c("blue", "turquoise", "midnightblue")
modNames <- substring(names(MEs), 3)  # removes "ME" prefix

# Calculate module membership (MM) for each gene (correlation with MEs)
geneModuleMembership <- as.data.frame(cor(gxData_t, MEs, use = "pairwise.complete.obs"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# Rename columns with "MM.<module>"
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Calculate gene significance (GS) for etiology
geneTraitSignificance <- as.data.frame(cor(gxData_t, metadata$etiology, use = "pairwise.complete.obs"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(metadata$etiology), sep = "")
names(GSPvalue) <- paste("p.GS.", names(metadata$etiology), sep = "")

# Plot
labels <- c("A)", "B)", "C)")
png("MM_vs_GS_top_modules.png", width = 2400, height = 1000, res = 300)
par(mfrow = c(1, length(modules)), mar = c(5, 5, 4, 2))

for (i in seq_along(modules)) {
  module <- modules[i]
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  
  # Base R scatter plot
  plot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = "Module Membership (MM)",
    ylab = "Gene Significance (GS) for DCM",
    main = paste0(toupper(substring(module,1,1)), substring(module,2)),
    col = module,
    cex = 1,
    cex.main = 1.3
  )
  
  # Add custom A)/B)/C) label in top-left, above main title
  title(
    main = labels[i],
    adj = 0,           # left-align
    col.main = "grey50",
    font.main = 2,
    cex.main = 1.5,
    line = 1.4
  )
}

dev.off()

# Build gene information table including only top modules --------------------#
geneInfo0 <- data.frame(
  Gene.ID = colnames(gxData_t),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)

# Find indices of top modules in modNames (to subset MM columns accordingly)
topModuleIndices <- match(modules, modNames)

# Add MM and p.MM columns for top modules only, keeping order consistent
for (idx in topModuleIndices) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, 
    geneModuleMembership[, idx], 
    MMPvalue[, idx]
  )
  names(geneInfo0) <- c(
    oldNames, 
    paste("MM.", modNames[idx], sep = ""),
    paste("p.MM.", modNames[idx], sep = "")
  )
}

# Filter geneInfo0 to keep only genes assigned to the top modules
geneInfo <- geneInfo0[geneInfo0$moduleColor %in% modules, ]

# Order genes by module color (following order in 'modules') and then by gene significance
geneOrder <- order(factor(geneInfo$moduleColor, levels = modules), -abs(geneInfo[[paste0("GS.", names(metadata$etiology))]]))
geneInfo <- geneInfo[geneOrder, ]

# Save to CSV file
write.csv(geneInfo, file = "geneInfo_topModules.csv", row.names = FALSE)

geneInfo.top50 <- geneInfo[1:50, c("Gene.ID", "GS.")]
write.csv(geneInfo.top50, file = "geneInfo_top50.csv", row.names = FALSE)

#-----------------------------------------------------------------------------#
# Cytoscape 
#-----------------------------------------------------------------------------#
RCy3::cytoscapePing()

if (!grepl("Installed", RCy3::getAppStatus("stringApp"))) {
  RCy3::installApp("stringApp")
}

# STRING network (BLUE MODULE) -----------------------------------------------#
genes.blue <- geneInfo[geneInfo$moduleColor == "blue", ]$Gene.ID

query <- format_csv(
  as.data.frame(genes.blue),
  col_names = FALSE,
  quote_escape = "double",
  eol = ",")

commandsPOST(
  paste0(
    'string protein query cutoff=0.4 ',
    'newNetName="Blue Module PPI Network" ',
    'query="', query, '" ',
    'limit=0'))

dataset <- read_excel("Data/DCM_diffexp_corr.xlsx")

loadTableData(
  dataset,
  data.key.column  = "Ensembl_GeneID",
  table.key.column = "query term")

# Network topology analysis & metrics ----------------------------------------#
analyzeNetwork(directed = FALSE)
node.table <- getTableColumns("node")

network.metrics <- node.table[, c(
  "display name",
  "query term",
  "Degree",
  "BetweennessCentrality",
  "ClusteringCoefficient"
)]

network.metrics <- network.metrics[
  order(network.metrics$Degree, decreasing = TRUE),]

# Visual Style set up --------------------------------------------------------#
RCy3::copyVisualStyle("default", "network_analysis_style")

# Node size = Degree (hub genes larger)
setNodeSizeMapping(
  table.column = "Degree",
  table.column.values = range(node.table$Degree, na.rm = TRUE),
  sizes = c(20, 120),
  mapping.type = "continuous",
  style.name = "network_analysis_style"
)

RCy3::setNodeColorMapping(
  table.column = "logFC",
  mapping.type = "continuous",
  colors = c("#67A9CF", "#FFFFFF", "#EF8A62"),
  table.column.values = c(
    min(dataset$logFC, na.rm = TRUE),
    0,
    max(dataset$logFC, na.rm = TRUE)
  ),
  style.name = "network_analysis_style"
)

# Node labels
RCy3::setNodeLabelMapping(
  table.column = "display name",
  style.name = "network_analysis_style"
)

RCy3::setVisualStyle("network_analysis_style")

# Hub Gene inspection ---------------------------------------------------------#
hist(network.metrics$Degree, breaks = 30, main = "Degree Distribution")
head(network.metrics, 10)

betweenness.genes <- network.metrics[order(network.metrics$BetweennessCentrality, decreasing = TRUE),]
head(betweenness.genes, 10)

# Network statistics ---------------------------------------------------------#
cat("\nNetwork Statistics:\n")
cat("  - Total nodes:", nrow(node.table), "\n")
cat("  - Average degree:", round(mean(node.table$Degree, na.rm = TRUE), 2), "\n")
cat("  - Network density:",
  round(mean(node.table$Degree, na.rm = TRUE) / (nrow(node.table) - 1), 3),
  "\n")
cat("  - Average clustering coefficient:",
  round(mean(node.table$ClusteringCoefficient, na.rm = TRUE), 3),
  "\n")

# Save outputs ---------------------------------------------------------------#
write.csv(
  network.metrics,
  file = "blue_module_network_metrics.csv",
  row.names = FALSE
)

exportImage(                              # image of Cytoscape network
  filename = "blue_module_network.png",
  type = "PNG",
  zoom = 200
)



#-----------------------------------------------------------------------------#
# GO Enrichment of WGCNA Genes
#-----------------------------------------------------------------------------#

# Convert Ensembl IDs to Entrez Gene IDs for GO analysis
genes.blue.entrez <- bitr(
    genes.blue,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )


cat("Genes mapped for GO enrichment:", nrow(genes.blue.entrez), "/", length(genes.blue), "\n\n")

go.bp <- enrichGO(
  gene = genes.blue.entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Convert to data frame for inspection
go.bp.res <- as.data.frame(go.bp)

cat("GO enrichment analysis results:\n") # Summary statistics
cat("  - Significant processes:", nrow(go.bp.res), "\n")

head(go.bp.res)

go.bp.res$Result

# Tree plot showing enriched pathways
num.bp <- min(20, nrow(go.bp.res))  # Do I need this line as well?

p <- treeplot(go.bp.sim, showCategory = 15) +
  theme(
    axis.text.y = element_blank(),   # <-- removes ONLY node labels
    axis.ticks.y = element_blank()   # optional, cleaner
  )

png("GO_BP_tree.png", width = 4000, height = 4000, res = 300)
print(p)
dev.off()


# Extract gene list ----------------------------------------------------------#
go_res <- as.data.frame(go.bp)
sig_go_terms <- go_res[go_res$p.adjust < 0.05, ]

# Extract genes from all significant GO terms
genes_entrez_list <- unlist(
  strsplit(sig_go_terms$geneID, split = "/"))

genes_entrez_unique <- unique(genes_entrez_list) # Remove duplicates


genes_entrez_unique <- genes_entrez_unique[genes_entrez_unique != "" & !is.na(genes_entrez_unique)]
genes_entrez_unique <- as.character(genes_entrez_unique)


# Convert Entrez IDs back to gene symbols for readability
sig_genes_symbols <- bitr(
  genes_entrez_unique,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)


# View significant genes linked to your enriched GO terms
print(sig_genes_symbols)