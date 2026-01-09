#=============================================================================#
# SVA.R                                                                       #
#	                                                                            #
# Surrogate Variable Analysis to correct for batch effects.                   #
#                                                                             #
# The SVA performed in the DiffExpressionAnalysis.R file of the MAGNet        #
# GitHub repository was used as a guideline for the analysis below.           #
# https://github.com/mpmorley/MAGNet/blob/master/bin/DiffExpressionAnalysis.R #
#                                                                             #
# Any code taken directly from the repository is specified.                   #
#                                                                             #
#=============================================================================#

# Load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(limma)
library(sva)
library(pcaMethods)

# Unsuccessful attempt at SVA batch correction :(

# ----------------------------------------------------------------------------#
# Downloaded already corrected data from GitHub repo previously mentioned
# Stored in workspace as R object (is RMD file)
# Uncomment line below  to write it as a CSV again

# write.csv(CPMS_SVA_corrected, "Data/CPMS_SVA_corrected.csv")
# ----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

gxData <- read.csv("Data/gxData_female.csv", as.is = T, row.names = 1)
metadata <- read.csv("Data/metadata_female.csv", as.is = T, row.names = 1)

#-----------------------------------------------------------------------------#
# PCA
# (Separate PCA to explore batch correction)
#-----------------------------------------------------------------------------#

pca_res <- pca(t(gxData), nPcs = 10, method = "nipals")

# For plotting
plot_pca <- cbind(data.frame(pca_res@scores), metadata)

pca1v2 <- ggplot(plot_pca, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = etiology, col=Library.Pool)) +
  theme_minimal() 

pca1v2

#-----------------------------------------------------------------------------#
# SVA
#-----------------------------------------------------------------------------#

# Check gxData and metadata samples are in the same order
all(colnames(gxData) == rownames(metadata))

# Create a Model with the variable of interest (etiology, whose variation we want to explain in our analysis)~
# Added covariates race and age
model <- model.matrix(~race+age+etiology, data = metadata)

# Create a Null model (baseline) without a variable of interest
model0 <- model.matrix(~1, data = metadata)  # Explains overall avg expression

# Taken from DiffExpressionAnalysis.R in repository, lines 20-30 --------------
# Variable names modified for clarity
#
# @Description
# Helper function to execute SVA
#
# data - count data (CPM)
# mod - model matrix
# null_mod - null model matrix
# n.sv - number of surrogate variables

svaBatchCor <- function(data, mod, null_mod, n.sv=NULL){
  # coerce data to matrix and transpose
  data <- as.matrix(data)
  Y <- t(data)
  
  # if n.sv not specified, estimate with leek method
  if(is.null(n.sv)){
    n.sv <- num.sv(data, mod, method="leek")
  }   
  
  # estimate surrogate variables
  svaRes <- svaseq(data, mod, null_mod, n.sv=n.sv)
  W <- svaRes$sv # this matrix is samples x SVs
  
  # matrix ordinary least squares regression
  # finds how much each gene's expression explained by random variation
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  
  # random variation subtracted 
  svaRes$corrected <- t(Y - W %*% alpha)
  
  return(svaRes)
}
# -----------------------------------------------------------------------------
sv <- svaBatchCor(gxData, model, model0, n.sv=24)

# Check which surrogate variables explain the most variance
sv.vars <- apply(sv.obj$sv, 2, var)
sv.rank <- order(sv.vars, decreasing = TRUE)

vars.df <- data.frame(
  SV = factor(paste0("SV", 1:ncol(sv.obj$sv)), levels = paste0("SV", 1:ncol(sv.obj$sv))),
  Variance = sv.vars)

ggplot(vars.df, aes(x = SV, y = Variance)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Variance Explained by Surrogate Variables", y = "Variance", x = "Surrogate Variable")

# Checking if the variances are indeed all the same (as seen in the scree plot above)
options(digits=15)
apply(sv.obj$sv, 2, var)
