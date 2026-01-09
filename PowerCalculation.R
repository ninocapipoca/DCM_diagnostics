#=============================================================================#
# PowerCalculation.R                                                          #
#	                                                                            #
#                                                                             #
# Start date: 09/01/2026      									                              #
#=============================================================================#

# Load packages ------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(gt) 
library(gtsummary)
library(pcaMethods)
library(limma)
library(PROPER)
library(Biobase)

#-----------------------------------------------------------------------------#
# RAW data import & Filtering
#-----------------------------------------------------------------------------#
# Set working directory according to the location of this script
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Import RAW data
RAWgxData_all <- read.csv("RawCounts.csv", as.is = T, row.names = 1)

# Filter RAWgxData_all -> keep only female DCM + NF samples
# NOTE - 'metadata' is previously created (DataPreparationExplorationAnalysis.R !)
RAWgxData <- RAWgxData_all[, colnames(RAWgxData_all) %in% rownames(metadata)]

# Remove lowly expressed genes
expr_samples <- rowSums(RAWgxData > 0)  # Number of samples each gene is expressed in
expr_cutoff <- 0.8 * ncol(RAWgxData)    # Create cutoff (80% of samples)
RAWgxData <- RAWgxData[expr_samples >= expr_cutoff, ]  # Filter RAWgxData with the cutoff

# Just making sure rownames(metadata) match colnames(RAWgxData)
all(colnames(RAWgxData) == rownames(metadata))  # Should be TRUE 

# Convert RAWgxData into an ExpressionSet (expected by estParam() later)
RAWgxData.mat <- as.matrix(RAWgxData)
exprSet <- ExpressionSet( 
  assayData = RAWgxData.mat,
  phenoData = AnnotatedDataFrame(metadata))

#-----------------------------------------------------------------------------#
# Setting Parameters
#-----------------------------------------------------------------------------#
# NOTE - The PROPER package simulates *new* data according to our study design; doesn't use our actual data
total_genes <- nrow(RAWgxData) # Total n genes
n1 <- 66                       # n DCM samples
n2 <- 89                       # n NF samples
logfc <- 0.5                   # Target log2 fold change (can change)
fdr_cut <- 0.5                 # FDR cutoff (can change)
DE_prop <- 0.05                # Expected proportion of DE genes
sim_num <- 100                 # n of simulations (can change: higher = longer, but more precise)

# Estimate simulation parameters from your raw count data
paramEst <- estParam(exprSet)

# Create a simulation object (for 2 groups, DCM + NF)
sim.obj <- RNAseq.SimOptions.2grp(
  ngenes = total_genes,
  seqDepth = paramEst$seqDepth,     # Sequencing depth
  lBaselineExpr = paramEst$lmeans,  # Typical expression level (log) in our samples
  lOD = paramEst$lOD,               # Overdispersion (log) per gene (RNAseq data tends to be)
  p.DE = DE_prop,
  lfc = logfc)

#-----------------------------------------------------------------------------#
# Wrapping Nina's DGE pipeline
#-----------------------------------------------------------------------------#
Nina_DGEpipe <- function(counts, n1, n2, fdr_cut) {
  
  # Simulated metadata
  metadata <- data.frame(
    etiology = factor(c(rep("DCM", n1), rep("NF", n2)))
  )
  
  design <- model.matrix(~ 0 + etiology, data = metadata)
  fit <- lmFit(counts, design)
  
  cont.matrix <- makeContrasts(
    DCMvsControl = etiologyDCM - etiologyNF,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  ebFit <- eBayes(fit2, trend = TRUE)
  
  dge_res <- topTable(
    ebFit,
    coef = "DCMvsControl",
    number = nrow(counts),
    sort.by = "none"
  )
  
  sig <- dge_res$adj.P.Val < fdr_cut
  return(sig)
}

#-----------------------------------------------------------------------------#
# Run simulations and calculate power
#-----------------------------------------------------------------------------#
power_vec <- numeric(sim_num) # vector to be filled with simulation results

for (i in seq_len(sim_num)) {
  sim <- simRNAseq(sim.obj, n1 = n1, n2 = n2)
  
  sig_genes <- Nina_DGEpipe(
    counts = sim$counts,
    n1 = n1,
    n2 = n2,
    fdr_cut = fdr_cut)
  
  true_DE <- sim$DEid
  
  power_vec[i] <- mean(sig_genes[true_DE], na.rm = TRUE)
}

mean(power_vec)
sd(power_vec)
quantile(power_vec, c(0.025, 0.975)) # 95% Monte Carlo uncertainty interval for power
