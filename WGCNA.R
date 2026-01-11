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
