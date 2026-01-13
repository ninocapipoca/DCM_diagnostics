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
# Pre-processing
#-----------------------------------------------------------------------------

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# Filter lowly-expressed genes (doesn't seem to have an effect)
keep <- rowSums(gxData > 1) >= 0.2 * ncol(gxData)
gxData_filtered <- gxData[keep, ]

# Sanity check -  zeroes and NaNs
sum(gxData == 0)
sum(is.na(gxData))



#-----------------------------------------------------------------------------
#
# Machine Learning
#
#-----------------------------------------------------------------------------

# Split dataset into training and test
set.seed(42)
n <- ncol(gxData)

gxData_t <- t(gxData)
resp_vect <- metadata$etiology[match(colnames(gxData), rownames(metadata))]
resp_vect <- factor(resp_vect)

# 70 train / 30 test split
index <- sample(1:n, .7*n)

x_train <- gxData_t[index,]
x_test <- gxData_t[-index,]

y_train <- resp_vect[index]
y_test <- resp_vect[-index]

#-----------------------------------------------------------------------------
# LASSO regression
#-----------------------------------------------------------------------------

# determine optimal values for lambda
# by default, 10-fold cross validation
lasso.fit <- cv.glmnet(x_train, y_train, 
                       alpha = 1, 
                       family = "binomial", # for classification
                       standardize = TRUE)

# Extract nonzero coefficients (important genes)
coeffs <- coef(lasso.fit, s = "lambda.1se")
key_genes <- rownames(coeffs)[coeffs[,1] != 0]
key_genes <- key_genes[key_genes != "(Intercept)"]

# Write to text file
cat(key_genes, file = "Data/key_genes_lasso.txt")

# Test set predictions
lasso.pred <- predict(lasso.fit,
                      newx = x_test,
                      s = "lambda.1se",
                      type = "response")

# Threshold 0.5
pred_class <- ifelse(lasso.pred > 0.5, "DCM", "NF")
pred_class <- factor(pred_class, levels = levels(y_test))

# Confusion matrix
table(Predicted = pred_class, Actual = y_test)





