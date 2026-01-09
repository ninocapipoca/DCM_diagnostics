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

#-----------------------------------------------------------------------------#
# Data import & parameter set-up
#-----------------------------------------------------------------------------#
# Set working directory according to the location of this script
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# NOTE - Must have loaded the gxData and metadata files previously.

# Setting parameters
# NOTE - The PROPER package simulates *new* data according to our study design, not use our actual data)
total_genes <- nrow(gxData) # Total number of genes
n1 <- 66                    # Number of DCM samples
n2 <- 89                    # Number of NF samples
logfc <- 0.5                # Target log2 fold change
fdr_cut <- 1.0              # FDR cutoff

args(estParam)
str(gxData)
class(gxData)
class(gxData.mat)

# Estimate simulation parameters from your raw count data
gxData.mat <- data.matrix(gxData)
class(gxData.mat) <- "matrix"
paramEst <- estParam(gxData.mat)

?simRNAseq # (Just checking functions!)
args(RNAseq.SimOptions.2grp)

# Create a simulation object (for 2  groups, DCM + NF) to feed to the simRNAseq() function later
sim.obj <- RNAseq.SimOptions.2grp(ngenes = total_genes,
                                  seqDepth = paramEst$seqDepth,
                                  nSamples = c(n1, n2),
                                  p.DE = 0.05,            # Expected proportion of DE genes
                                  lfc = logfc)

# Run RNA-seq power simulations using limma-voom
simResults <- simRNAseq(simOpts, testMethod = "limma_voom", nsim = sim_num, fdr = fdr_cut)

# Summarize power
summaryPower(simResults)

# Optional: plot power
plotPower(simResults)

