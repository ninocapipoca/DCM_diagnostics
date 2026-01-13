library(glmnet)
library(dplyr)

# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load data
data <- read.csv("Data/CPMS_SVA_corrected.csv", row.names = 1, header = TRUE)
metadata_SVA <- read.csv("Data/phenoData_SVA.csv", row.names = 1, header = TRUE)

metadata <- metadata_SVA |> select(1:16)

# Split dataset into training and test
set.seed(42)
n <- ncol(data)

train <- sample(1:n, size = round(0.7 * num_samples)) # use 70% as trainining set
#test <- 

