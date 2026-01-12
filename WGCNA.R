#=============================================================================#
# WGCNA.R                                                                     #
#	                                                                            #
# Start date: 11/01/2026      									                              #
#=============================================================================#

# Load packages ------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(impute)
library(preprocessCore)
library(WGCNA)
library(readr)
library(readxl)
library(RCy3)
library(rstudioapi)

allowWGCNAThreads()

#-----------------------------------------------------------------------------#
# Data import (same as in DataPreparationExploratoryAnalysis.R)
#-----------------------------------------------------------------------------#
# Set working directory according to the location of this script
setwd("/Users/mikiverme/Desktop/DCM_diagnostics")

# Load gene expression data and participant information (metadata)
# NOTE - Gene expression dataset contains log2-transformed CPM
gxData_all <- read.table("MAGNET_GX_2025/MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T, row.names = 1, header = TRUE)
metadata_all <- read.csv("MAGNET_GX_2025/MAGNET_SampleData_18112022.csv", as.is = T, row.names = 1)

# Filter out only female participants who either have DCM or are healthy
metadata <- metadata_all |> 
  filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))

# Filter gxData_all to include only participants also in metadata
gxData <- gxData_all[, colnames(gxData_all) %in% rownames(metadata)]

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# TODO - Possibly remove outliers (after talking to Aaron on Monday)

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
metadata <- mutate_all(metadata, function(x) as.numeric(as.character(x)))

# Display trait summary
summary(metadata)

# Clean up memory
collectGarbage()

# Sample clustering and trait visualisation ----------------------------------#

# Hierarchical clustering of samples based on gene expression
sampleTree <- hclust(dist(gxData_t), method = "average")

# Convert traits to color representation for visualization
# White = low value, red = high value, grey = missing
traitColors <- numbers2colors(metadata, signed = FALSE) # NOTE - unsigned!!

# Plot and save dendrogram with trait heatmap
dendrogram_path <- file.path("/Users/mikiverme/Desktop/DCM_diagnostics", "dendrogram_heatmap.png")
png(filename = dendrogram_path, width = 1200, height = 900, res = 150)

plotDendroAndColors(
  sampleTree, 
  traitColors,
  groupLabels = names(metadata), 
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


# Visualise power selection
# TODO - Fix figures ! 
par(mfrow = c(1, 2))

# Plot 1: Scale-free topology fit (R²) vs. power
plot(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (signed R²)",
  type = "n", 
  main = "Scale Independence"
)

text(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
  labels = powers,
  cex = 0.9,
  col = "red"
)

abline(h = 0.85, col = "red", lty = 2) # Horizontal line at R² = 0.85 (recommended threshold)


# Plot 2: Mean connectivity vs. power
plot(
  sft$fitIndices[, 1], 
  sft$fitIndices[, 5], 
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity", 
  type = "n", 
  main = "Mean Connectivity"
)
text(
  sft$fitIndices[, 1], 
  sft$fitIndices[, 5], 
  labels = powers, 
  cex = 0.9,
  col = "red"
)


# Build co expression network ------------------------------------------------#
# Construct the network and identify modules
net <- blockwiseModules(
  gxData_t, 
  power = 6,                      # Soft-thresholding power
  TOMType = "unsigned",           # Topological Overlap Matrix type
  minModuleSize = 30,             # Minimum module size
  reassignThreshold = 0,          # Gene reassignment threshold
  mergeCutHeight = 0.25,          # Module merging threshold
  numericLabels = TRUE,           # Use numbers instead of colors initially
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "expTOM", 
  verbose = 3
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

# GENE BLOCK 2
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

# GENE BLOCK 3
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


# TODO - Save the dendrograms


# Module Trait Relationships -------------------------------------------------#

# Calculate module eigengenes (MEs) - first principal component of each module
nGenes <- ncol(gxData_t)
nSamples <- nrow(gxData_t)

MEs0 <- moduleEigengenes(gxData_t, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate MEs with clinical traits
moduleTraitCor <- cor(MEs, metadata, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

round(moduleTraitCor, 2) # Display correlation matrix

# Create text matrix with correlations and p-values
textMatrix <- paste(
  signif(moduleTraitCor, 2), 
  "\n(",
  signif(moduleTraitPvalue, 1), 
  ")", 
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)

# Set margins for plot
par(mar = c(8, 8.5, 3, 3))

# Create heatmap
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(metadata),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  cex.lab.y = 0.6,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)

# TODO - Identify what variables (from metadata) do we want to keep in our study
#        (it defines how many columns in the heatmap), label the variables in the
#        heatmap, and save the heatmap


# Gene module membership and significance ------------------------------------#

# Define etiology as the trait of interest
eti <- as.data.frame(metadata$etiology)
names(eti)
modNames <- substring(names(MEs), 3)

# Calculate module membership (MM) for each gene
geneModuleMembership <- as.data.frame(cor(gxData_t, MEs, use = "pairwise.complete.obs"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Calculate gene significance (GS) for etiology
geneTraitSignificance <- as.data.frame(cor(gxData_t, metadata$etiology, use = "pairwise.complete.obs"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(metadata$etiology), sep = "")
names(GSPvalue) <- paste("p.GS.", names(metadata$etiology), sep = "")

# Visualise the MM vs GS
modules <- c("brown", "black", "lightyellow")
par(mfrow = c(1, 3))

for (module in modules) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership (MM) in", module, "module"),
    ylab = "Gene Significance (GS) for disease",
    main = paste("MM vs. GS\n"),
    cex.main = 1.2, 
    cex.lab = 1, 
    cex.axis = 1, 
    col = module
  )
}

# TODO - fix formatting of the graphs! and know how to interpret them

# Create comprehensive gene information table
geneInfo0 <- data.frame(
  Gene.ID = colnames(gxData_t),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)

# Order modules by their significance for disease
modOrder <- order(-abs(cor(MEs, metadata$etiology, use = "pairwise.complete.obs")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, 
    geneModuleMembership[, modOrder[mod]], 
    MMPvalue[, modOrder[mod]]
  )
  names(geneInfo0) <- c(
    oldNames, 
    paste("MM.", modNames[modOrder[mod]], sep = ""),
    paste("p.MM.", modNames[modOrder[mod]], sep = "")
  )
}

# Order genes by module color, then by gene significance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.))
geneInfo <- geneInfo0[geneOrder, ]

# Save to CSV file
write.csv(geneInfo, file = "geneInfo.csv", row.names = FALSE)

#-----------------------------------------------------------------------------#
# Cytoscape 
#-----------------------------------------------------------------------------#

# Test connection to Cytoscape and install stringApp if needed ---------------#

RCy3::cytoscapePing() # Test connection to Cytoscape

if (!grepl("status: Installed", RCy3::getAppStatus("stringApp"))) {
  RCy3::installApp("stringApp")
  cat("✓ stringApp installed\n")
}

# Create PPI Network ---------------------------------------------------------#

# BLACK MODULE gene selection
genes.black <- geneInfo[geneInfo$moduleColor == "black", ]$Gene.ID

query <- format_csv(          # Format gene list for STRING query
  as.data.frame(genes.black), 
  col_names = FALSE, 
  quote_escape = "double", 
  eol = ",")

commandsPOST(                 # Query STRING database through Cytoscape
  paste0(
    'string protein query cutoff=0.4 newNetName="Black Module PPI Network" query="',
    query,
    '" limit=0'))

dataset <- read_excel("DCM_diffexp_copy.xlsx")    # Load differential expression data onto network
loadTableData(
  dataset, 
  data.key.column = "Ensembl_GeneID", 
  table.key.column = "query term")

# BLACK MODULE Adding expression heatmap -------------------------------------#

RCy3::copyVisualStyle("default", "my_style_heatmap2")

# Add heatmap showing male and female differential expression
RCy3::setNodeCustomHeatMapChart("logFC",
  slot = 2, 
  style.name = "my_style_heatmap2", 
  colors = c("#EF8A62", "#FFFFFF", "#67A9CF"),  # Red-White-Blue
  range = c(min(dataset$logFC), max(dataset$logFC)))

RCy3::setNodeLabelMapping("display name", style.name="my_style_heatmap2")

# Apply visual style
RCy3::setVisualStyle("my_style_heatmap2")

# BLACK MODULE Analyze network topology --------------------------------------#

# This calculates degree, betweenness centrality, clustering coefficient, etc.
analyzeNetwork(directed = FALSE)

# Get node table with all network properties
node.table <- getTableColumns(table = "node")

# Extract key network metrics
network.metrics <- node.table[, c(
  "display name",
  "Degree",
  "BetweennessCentrality", 
  "ClusteringCoefficient"
)]

# Sort by degree to find hub genes
network.metrics <- network.metrics[order(network.metrics$Degree, decreasing = TRUE), ]
hist(network.metrics$Degree, breaks=30, main = "Degree distribution")

# Display top 10 hub genes
print(head(network.metrics, 10))

# Get network statistics
cat("\nNetwork Statistics:\n")
cat("  - Total nodes:", nrow(node.table), "\n")
cat("  - Average degree:", round(mean(node.table$Degree, na.rm = TRUE), 2), "\n")
cat("  - Network density:", round(mean(node.table$Degree, na.rm = TRUE) / (nrow(node.table) - 1), 3), "\n")
cat("  - Average clustering coefficient:", 
    round(mean(node.table$ClusteringCoefficient, na.rm = TRUE), 3), "\n")

# Save network metrics to file
write.csv(network.metrics, file = "black_module_network_metrics.csv", row.names = FALSE)

# Create custom visual style
RCy3::copyVisualStyle("default", "network_analysis_style")

# Map node size to Degree (larger nodes = more connections)
setNodeSizeMapping(
  table.column = "Degree",
  table.column.values = c(min(node.table$Degree), max(node.table$Degree)),
  sizes = c(30, 150),
  mapping.type = "continuous",
  style.name = "network_analysis_style"
)

# Add expression heatmap
RCy3::setNodeCustomHeatMapChart(
  c("DCMvsControl_male_logFC", "DCMvsControl_female_logFC"), 
  slot = 2, 
  style.name = "network_analysis_style", 
  colors = c("#EF8A62", "#FFFFFF", "#67A9CF"),
  range = c(-1, 1)
)

RCy3::setNodeLabelMapping("display name", style.name="network_analysis_style")


# Apply visual style
RCy3::setVisualStyle("network_analysis_style")

cat("✓ Visual style applied!\n")
cat("  - Node size = Degree\n")
cat("  - Node color = Betweenness Centrality\n")
cat("  - Heatmap = Male/Female expression changes\n")

# Export network image
exportImage(
  filename = "black_module_network.png",
  type = "PNG",
  zoom = 200
)


