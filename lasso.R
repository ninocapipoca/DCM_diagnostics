library(glmnet)
library(dplyr)

# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load data
data <- read.csv("Data/CPMS_SVA_corrected.csv", row.names = 1, header = TRUE)
metadata_SVA <- read.csv("Data/phenoData_SVA.csv", row.names = 1, header = TRUE)

#-----------------------------------------------------------------------------
# Filter dataset
#-----------------------------------------------------------------------------

# Filter out SVA columns
metadata <- metadata_SVA |> select(1:16)

# Filter out only female participants who either have DCM or are healthy
metadata <- metadata |> 
  filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))

# Filter gene data to include only participants also in metadata
gxData <- data[, colnames(data) %in% rownames(metadata)]

#-----------------------------------------------------------------------------
# Preprocessing
#-----------------------------------------------------------------------------

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# Filter lowly-expressed genes
keep <- rowSums(gxData > 5) >= 0.2 * ncol(gxData)
gxData_filtered <- expr_data[keep_genes, ]


#-----------------------------------------------------------------------------
#
# Machine Learning
#
#-----------------------------------------------------------------------------

# Split dataset into training and test
set.seed(42)
n <- ncol(gxData)

gxData_t <- t(gxData)
resp_vect <- colnames(gxData_t)

# 70 train / 30 test split
index <- sample(1:n, .7*n)

x_train <- gxData_t[index,]
x_test <- gxData_t[-index]

y_train <- resp_vect[index]
y_test <- resp_vect[-index]

#-----------------------------------------------------------------------------
# LASSO regression
#-----------------------------------------------------------------------------

# determine optimal values for lambda
lasso.fit <- cv.glmnet(x_train, y_train, alpha=1, standardize=TRUE)




