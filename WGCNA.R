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
library(missMethods)

allowWGCNAThreads()


#-----------------------------------------------------------------------------#
# Data import (same as in DataPreparationExploratoryAnalysis.R)
#-----------------------------------------------------------------------------#
# Set working directory according to the location of this script
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load gene expression data and participant information (metadata)
# NOTE - *SVA-corrected* gene expression dataset (log2-transformed CPM)
gxData_all <- read.csv("Data/CPMS_SVA_corrected.csv", as.is = T, row.names = 1) |> 
  log2(gxData_all + 1) |>  #Log2 transform
  impute_median(gxData_all, type = "columnwise")
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

meta.traits <- metadata[, c("etiology", "age", "race", "LVEF", "VTVF", "Hypertension", "Diabetes")]
traitColors <- numbers2colors(meta.traits, signed = FALSE) # Convert those traits to colors

# Plot and save dendrogram with trait heatmap
dendrogram_path <- file.path("/Users/mikiverme/Desktop/DCM_diagnostics", "dendrogram_heatmap.png")
png(filename = dendrogram_path, width = 1200, height = 900, res = 150)

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

targetPower.test <- 4
meanConnectivity <- sft$fitIndices[               # Test of mean connectivity for power '4'
  sft$fitIndices[, "Power"] == targetPower.test,
  "mean.k."
]

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
  power = 4,                      # Soft-thresholding power
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


# Module Trait Relationships ----------------------CORRECTED-------------------#

# Calculate module eigengenes (MEs) - first principal component of each module
nGenes <- ncol(gxData_t)
nSamples <- nrow(gxData_t)

MEs0 <- moduleEigengenes(gxData_t, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Select traits of interest explicitly
traits_of_interest <- metadata[, c("etiology", "age", "race", "LVEF", "Diabetes", "Hypertension")]  # adjust as needed

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

# Subset the top 6 ranked modules (for later)
top5.modules <- ranked.modules$Module[c(1, 2, 4, 5, 6)]


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
yLabels_sub <- top5.modules
ySymbols_sub <- top5.modules

png("top5_modules_heatmap.png", width = 1000, height = 800, res = 150)
par(mar = c(8, 8.5, 3, 3))  # Adjust margins for your plot

labeledHeatmap(
  Matrix = moduleTraitCor_sub,
  xLabels = colnames(moduleTraitCor_sub),
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

# TODO - Label the correlations legend (right, red/white/blue)

# Gene module membership and significance ------------------------------------#

# Define top modules (from  ranked.modules output)
modules <- c("blue", "turquoise", "darkslateblue")
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

# Plot MM vs GS only for the top modules
par(mfrow = c(1, length(modules)))

for (module in modules) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership (MM) in", module, "module"),
    ylab = "Gene Significance (GS) for disease",
    main = paste("MM vs. GS for", module),
    cex.main = 1.2, 
    cex.lab = 1, 
    cex.axis = 1, 
    col = module
  )
}

# Build gene information table including only top modules
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


#-----------------------------------------------------------------------------#
# Cytoscape 
#-----------------------------------------------------------------------------#

# Test connection to Cytoscape and install stringApp if needed

RCy3::cytoscapePing() # Test connection to Cytoscape

if (!grepl("status: Installed", RCy3::getAppStatus("stringApp"))) {
  RCy3::installApp("stringApp")
  cat("✓ stringApp installed\n")
}

# Create PPI Network ---------------------------------------------------------#

# BLUE MODULE gene selection
genes.blue <- geneInfo[geneInfo$moduleColor == "blue", ]$Gene.ID

# Format gene list for STRING database query via Cytoscape
query <- format_csv(
  as.data.frame(genes.blue), 
  col_names = FALSE, 
  quote_escape = "double", 
  eol = ","
)

# Query STRING protein-protein interaction network for blue module genes
commandsPOST(
  paste0(
    'string protein query cutoff=0.4 newNetName="Blue Module PPI Network" query="',
    query,
    '" limit=0'
  )
)

# Load differential expression dataset (female DCM vs female control)
dataset <- read_excel("Data/DCM_diffexp_corr.xlsx")

# Load dataset onto Cytoscape network nodes, matching by Ensembl Gene ID
loadTableData(
  dataset, 
  data.key.column = "Ensembl_GeneID", 
  table.key.column = "query term"
)

# Visual style setup ---------------------------------------------------------#

# Default Cytoscape style (base for customisation)
RCy3::copyVisualStyle("default", "network_analysis_style")

# Map node size to Degree centrality (nodes with more connections appear larger)
setNodeSizeMapping(
  table.column = "Degree",
  table.column.values = c(min(node.table$Degree), max(node.table$Degree)),
  sizes = c(30, 150),
  mapping.type = "continuous",
  style.name = "network_analysis_style"
)

# Map node color continuously to female log fold change (DCM vs control)
# Blue = downregulated, White = no change, Red = upregulated
RCy3::setNodeColorMapping(
  table.column = "logFC",
  mapping.type = "continuous",
  colors = c("#67A9CF", "#FFFFFF", "#EF8A62"),  # Blue - White - Red
  style.name = "network_analysis_style",
  table.column.values = c(
    min(dataset$logFC, na.rm=TRUE), 
    0, 
    max(dataset$logFC, na.rm=TRUE)
  )
)

# Map node labels to gene names for clarity
RCy3::setNodeLabelMapping("display name", style.name="network_analysis_style")

# Apply the customized visual style to the network
RCy3::setVisualStyle("network_analysis_style")

cat("✓ Visual style applied!\n")
cat("  - Node size = Degree (connectivity)\n")
cat("  - Node color = Female logFC (expression change)\n")

# --- Network topology analysis ---------------------------------------------

# Calculate network topology metrics: degree, betweenness centrality, clustering coefficient
# Retrieve node attribute table with calculated metrics
analyzeNetwork(directed = FALSE)
node.table <- getTableColumns(table = "node")

# Extract key metrics for downstream analysis
network.metrics <- node.table[, c(
  "display name",
  "Degree",
  "BetweennessCentrality", 
  "ClusteringCoefficient"
)]

# Sort nodes by degree to identify hub genes
network.metrics <- network.metrics[order(network.metrics$Degree, decreasing = TRUE), ]

# Visualize degree distribution across the network
hist(network.metrics$Degree, breaks=30, main = "Degree Distribution")

# Display top 10 hub genes based on degree
print(head(network.metrics, 10))

# Print overall network statistics
cat("\nNetwork Statistics:\n")
cat("  - Total nodes:", nrow(node.table), "\n")
cat("  - Average degree:", round(mean(node.table$Degree, na.rm = TRUE), 2), "\n")
cat("  - Network density:", round(mean(node.table$Degree, na.rm = TRUE) / (nrow(node.table) - 1), 3), "\n")
cat("  - Average clustering coefficient:", round(mean(node.table$ClusteringCoefficient, na.rm = TRUE), 3), "\n")

# Save network metrics to CSV file for record keeping
write.csv(network.metrics, file = "blue_module_network_metrics.csv", row.names = FALSE)

# --- Export network image ---------------------------------------------------

exportImage(
  filename = "blue_module_network.png",
  type = "PNG",
  zoom = 200
)
