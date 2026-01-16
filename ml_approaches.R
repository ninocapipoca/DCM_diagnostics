.packages <- c("dplyr", "pROC", "caret", "glmnet", "missMethods",
               "ggplot2")
lapply(.packages, require, character.only = TRUE)

# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

if (exists("CPMS_SVA_corrected")){
  # In case you download the .RMD file from GitHub and open it here
  data <- CPMS_SVA_corrected
} else {
  tryCatch({
    # Try loading gene expression data from corrected values
    data <- read.csv("Data/CPMS_SVA_corrected.csv", row.names = 1, header = TRUE)
  }, warning = function(w) {
    cat("A warning occurred:", conditionMessage(w), "\n")
    print("Add the file to Data or download from https://github.com/mpmorley/MAGNet")
  })
}

tryCatch({
  # Try loading metadata
  metadata_SVA <- read.csv("Data/phenoData_SVA.csv", row.names = 1, header = TRUE)
}, warning = function(w) {
  cat("A warning occurred:", conditionMessage(w), "\n")
  print("Try downloading the file & adding it to the Data folder: https://www.dropbox.com/s/eihem5fbnkg7bpm/phenoData.csv?dl=0")
})

#-----------------------------------------------------------------------------
# Preprocessing
#-----------------------------------------------------------------------------

# Log2 normalize
data <- log2(data + 1)

# Filter out SVA columns
metadata <- metadata_SVA |> dplyr::select(1:16)

.female <- TRUE # For exploring dataset further, if necessary
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
# Cleanup
#-----------------------------------------------------------------------------

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# Sanity check -  zeroes and NaNs
sum(gxData == 0)
sum(is.na(gxData))

# LASSO REGRESSION checks
# Looking at distributions of top & "meaningless" genes -------
# ENSG00000273270   
# ENSG00000273271

# ENSG00000163513, ENSG00000144895 
#gxData <- gxData[c("ENSG00000273270", "ENSG00000273271"),]

# Histogram of distribution
# ENSG00000163513
# .top <-  gxData["ENSG00000205795",]
# hist(as.numeric(.top))
# 
# .top_NF <- .top[, colnames(.top) %in% rownames(metadata[metadata$etiology == "NF",])]
# hist(as.numeric(.top_NF))
# 
# .top_DCM <- .top[, colnames(.top) %in% rownames(metadata[metadata$etiology == "DCM",])]
# hist(as.numeric(.top_DCM))
# ------------------------------------ 

#-----------------------------------------------------------------------------
#
# Machine Learning
#
#-----------------------------------------------------------------------------
# Transpose
gxData_t <- t(gxData)

# Split dataset into training and test
set.seed(42)
n <- nrow(gxData_t)

resp_vect <- metadata$etiology[match(rownames(gxData_t), rownames(metadata))]
resp_vect <- factor(resp_vect, levels = c("NF", "DCM"))

# 70 train / 30 test split
#index <- sample(1:n, .7*n) - for non stratified

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

# determine optimal values for lambda
# by default, 10-fold cross validation
lasso.fit <- cv.glmnet(x_train, y_train,
                       alpha = 1,
                       family = "binomial", # for classification
                       standardize = TRUE)

plot(lasso.fit)

# Extract nonzero coefficients (important genes)
coeffs <- coef(lasso.fit, s = "lambda.1se")
coeffs_nonzero <- coeffs[coeffs[,1] != 0]
key_genes <- rownames(coeffs)[coeffs[,1] != 0]
key_genes <- key_genes[key_genes != "(Intercept)"]

values <- coeffs_nonzero[2:length(coeffs_nonzero)]

key_genes_df <- data.frame(Gene = key_genes, lambda = values)
key_genes_df <- key_genes_df[order(abs(key_genes_df$lambda)), ]


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

# Check histogram & confusion matrix
summary(lasso.pred)
hist(lasso.pred, breaks = 20)

confusionMatrix(pred_class, y_test)

# Check ROC
roc_obj <- roc(y_test, as.numeric(lasso.pred))
auc(roc_obj)
plot(roc_obj)

###### WARNING - GENAI! ----------------------------------
# Used for permutation testing

# 1. Define the number of permutations
n_permutations <- 30
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

#-----------------------------------------------------------------------------
# Random forest
#-----------------------------------------------------------------------------

control <- trainControl(method = "cv", number = 5)

RF.model <- train(x = x_train,
                  y = y_train,
                  method = "rf",
                  trControl = control,
                  ntree = 100)

RF.model

# Evaluate performance on test set
SVM.predict <- predict(SVM.model, newdata = x_test)
confusionMatrix(SVM.predict, y_test)

# Check histogram
hist(as.numeric(SVM.predict), breaks = 20)

# Check ROC plot
roc_obj <- roc(y_test, as.numeric(SVM.predict))
auc(roc_obj)
plot(roc_obj)


#-----------------------------------------------------------------------------
# SVM
# There is something wrong here
#-----------------------------------------------------------------------------

# Specify info for training
# classProbs is class probabilities
control <- trainControl(method="cv",number=5, classProbs = TRUE)

SVM.model <- train(x = x_train,
                   y = y_train,
                   method = "svmLinear2", # linear kernel
                   tuneLength = 5, # granularity of tuning param grid
                   preProc = c("center","scale"), # preprocessing
                   metric='Accuracy', 
                   trControl=control)

SVM.model

# Evaluate performance on test set
SVM.predict <- predict(SVM.model, newdata = x_test)
confusionMatrix(SVM.predict, y_test)

# Check histogram
hist(as.numeric(SVM.predict), breaks = 20)

# Check ROC plot
roc_obj <- roc(y_test, as.numeric(SVM.predict))
auc(roc_obj)
plot(roc_obj)

###### WARNING - GENAI! ----------------------------------
# Used for permutation testing. Adjusted for SVM.

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
  fit_shuffled <- train(x = x_train,
                        y = y_train,
                        method = "svmLinear2",
                        tuneLength = 5,
                        preProc = c("center","scale"),
                        metric='Accuracy', 
                        trControl=control)
  
  # Predict on the (unshuffled) test set
  pred_shuffled <- predict(fit_shuffled, newdata = x_test)
  
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


