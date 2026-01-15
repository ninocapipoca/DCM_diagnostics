library(glmnet)
library(dplyr)
library(missMethods)
library(pROC)
library(caret)

# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load data
data <- read.csv("Data/CPMS_SVA_corrected.csv", row.names = 1, header = TRUE)
metadata_SVA <- read.csv("Data/phenoData_SVA.csv", row.names = 1, header = TRUE)

#-----------------------------------------------------------------------------
# Filter dataset
#-----------------------------------------------------------------------------

# Log2 normalize
data <- log2(data + 1)

# Filter out SVA columns
metadata <- metadata_SVA |> select(1:16)


.female <- TRUE # CHANGE THIS FLAG WHEN TESTING
if (.female){
  # Filter out only female participants who either have DCM or are healthy
  metadata <- metadata |> 
    filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))
} else {
  # Filter out everything that is not DCM or healthy
  metadata <- metadata |> 
    filter(etiology == "DCM" | etiology == "NF")
}


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

# Apply median imputation
x_train <- impute_median(x_train, type = "columnwise")
x_test <- impute_median(x_test, type = "columnwise")

sum(is.na(x_train))
sum(is.na(x_test))

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

###### WARNING - GENAI!!!!!!! 
# 1. Define the number of permutations
n_permutations <- 20
null_aucs <- numeric(n_permutations)

# 2. Get your "True" AUC first (from your original code)
# Assuming roc_obj is already calculated from your script
true_auc <- as.numeric(auc(roc_obj))

set.seed(123) # For reproducibility of the shuffles
for (i in 1:n_permutations) {
  
  # SHUFFLE the training labels
  y_train_shuffled <- sample(y_train)
  
  # Fit model on shuffled labels
  fit_shuffled <- cv.glmnet(x_train, y_train_shuffled, 
                            alpha = 1, 
                            family = "binomial", 
                            standardize = TRUE)
  
  # Predict on the (unshuffled) test set
  pred_shuffled <- predict(fit_shuffled, 
                           newx = x_test, 
                           s = "lambda.1se", 
                           type = "response")
  
  # Calculate AUC for this null run
  roc_shuffled <- roc(y_test, as.numeric(pred_shuffled), quiet = TRUE)
  null_aucs[i] <- as.numeric(auc(roc_shuffled))
  
  # Optional: Print progress
  if(i %% 10 == 0) cat("Completed", i, "permutations...\n")
}

# 3. Visualize the results
hist(null_aucs, main = "Permutation Test (Null Distribution of AUC)",
     xlab = "AUC with Shuffled Labels", xlim = c(0, 1), col = "lightblue")
abline(v = true_auc, col = "red", lwd = 2, lty = 2)
text(true_auc, 0.5, "True AUC", pos = 4, col = "red")

