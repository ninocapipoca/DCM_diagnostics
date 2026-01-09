# Load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(limma)
library(sva)


#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

gxData <- read.csv("Data/gxData_female.csv", as.is = T, row.names = 1)
metadata <- read.csv("Data/metadata_female.csv", as.is = T, row.names = 1)

#-----------------------------------------------------------------------------#
# SVA
#-----------------------------------------------------------------------------#
# NOTE - The results are very peculiar and I don't know if it's done correctly.. I will ask at next Q&A

# Check gxData and metadata samples are in the same order
all(colnames(gxData) == rownames(metadata))

# Create a Model with the variable of interest (etiology, whose variation we want to explain in our analysis)
model <- model.matrix(~ etiology, data = metadata)

# Create a Null model (baseline) without a variable of interest
model0 <- model.matrix(~ 1, data = metadata)  # Explains overall avg expression

# Create a surrogate variable (unknown factors; separates known variation from total variation and compares to baseline)
# Will be a matrix with rows = samples, cols = unknown factors (+ value per sample) that cause variability in the data
sv.obj <- sva(as.matrix(gxData), model, model0)
modelSV <- cbind(model, sv.obj$sv) # NOTE - the number of SVs is very high, so we check what's up in the next lines

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