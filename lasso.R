library(glmnet)
library(dplyr)
library(missMethods)
library(pROC)
library(caret)

# UNCORRECTED!!!!! DATA VERSION
# Using uncorrected data to see if SVA is source of leakage

# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load data
data <- read.csv("Data/MAGNET_GX_2025/gxData_female.csv", row.names = 1, header = TRUE)
metadata <- read.csv("Data/MAGNET_GX_2025/metadata_female.csv", row.names = 1, header = TRUE)

#-----------------------------------------------------------------------------
# Filter dataset
#-----------------------------------------------------------------------------

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# Sanity check -  zeroes and NaNs
sum(gxData == 0)
sum(is.na(gxData))

# Transpose
gxData_t <- t(gxData)



#-----------------------------------------------------------------------------
#
# Machine Learning
#
#-----------------------------------------------------------------------------

# Split dataset into training and test
set.seed(42)
n <- nrow(gxData_t)

resp_vect <- metadata$etiology[match(rownames(gxData_t), rownames(metadata))]
resp_vect <- factor(resp_vect, levels = c("NF", "DCM"))

# 70 train / 30 test split
#index <- sample(1:n, .7*n)
# stratified sampling to avoid one class being underrepresented
index <- createDataPartition(resp_vect, p = 0.7, list = FALSE) 

x_train <- gxData_t[index,]
x_test <- gxData_t[-index,]

y_train <- resp_vect[index]
y_test <- resp_vect[-index]

#-----------------------------------------------------------------------------
# LASSO regression
#-----------------------------------------------------------------------------

# # Another attempt at Lasso regression
# # This time with repeated CV 
# train_control <- trainControl(
#   method = "repeatedcv",
#   number = 10,             
#   repeats = 5,           
#   classProbs = TRUE,
#   summaryFunction = twoClassSummary
# )
# 
# lasso_model <- train(
#   x = x_train,
#   y = y_train,
#   method = "glmnet",
#   metric = "ROC",  # optimize for AUC
#   trControl = train_control,
#   tuneLength = 5
# )
# 
# lasso_model$results
# 
# pred_test <- predict(lasso_model, newdata = x_test)
# confusionMatrix(pred_test, y_test)


# determine optimal values for lambda
# by default, 10-fold cross validation
lasso.fit <- cv.glmnet(x_train, y_train,
                       alpha = 1,
                       family = "binomial", # for classification
                       standardize = TRUE)

plot(lasso.fit)

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

# Check histogram
summary(lasso.pred)
hist(lasso.pred, breaks = 20)

# Check ROC
roc_obj <- roc(y_test, as.numeric(lasso.pred))
auc(roc_obj)
plot(roc_obj)

