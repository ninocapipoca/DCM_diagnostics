#==============================================================================
# ml_approaches.R                                                             #
#																	                                       		  #
# Implementing various classifiers which provide rankings for genes,          #
# allowing us to obtain a list of most relevant genes in DCM.                 #
#=============================================================================#
# Flags
.female <- FALSE # Decide which part of the data to run on
.permute <- FALSE # Run permutation tests or not
.combine <- FALSE # Combine SVM & lasso plots, summary stats

# Choose which algorithms to run (default F, T, T)
.RF <- FALSE
.SVM <- TRUE
.lasso <- FALSE


# Load packages
.packages <- c("dplyr", "pROC", "caret", "glmnet", "missMethods",
               "ggplot2", "e1071", "gt", "gtsummary", "patchwork")
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

# ----------------------------------------------------------------------------
# Define ggplot theme
# ----------------------------------------------------------------------------

theme_pretty <- function() {
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.background = element_rect(fill = "#f0f0f0"),
    panel.grid.major = element_line(colour = "#c5c5c5"),
    panel.grid.minor = element_blank()
  )
}

#-----------------------------------------------------------------------------
# Preprocessing
#-----------------------------------------------------------------------------

# Log2 normalize
data <- log2(data + 1)

# Filter out SVA columns
metadata <- metadata_SVA |> dplyr::select(1:16)

if (.female){
  s <- "F"
  # Filter: only female participants who either have DCM or are healthy
  metadata <- metadata |> 
    filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))
} else {
  s <- "M"
  # Filter: only male participants
  metadata <- metadata |> 
    filter(gender == "Male" & (etiology == "DCM" | etiology == "NF"))
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
# LASSO classifier
#-----------------------------------------------------------------------------
if (.lasso){
  print("Running LASSO")
  # determine optimal values for lambda
  # by default, 10-fold cross validation
  lasso.fit <- cv.glmnet(x_train, y_train,
                         alpha = 1,
                         family = "binomial", # for classification
                         standardize = TRUE)
  
  plot(lasso.fit)
  
  # PRETTYPLOTS - Lasso performance and binomialdev
  lasso_cv_data <- data.frame(
    lambda = lasso.fit$lambda,
    cvm = lasso.fit$cvm,
    cvsd = lasso.fit$cvsd,
    num_features <- colSums(lasso.fit$glmnet.fit$beta != 0) # nonzero coeffs
  )
  
  lasso_perf <- ggplot(lasso_cv_data, aes(x = num_features, y = 1 - cvm)) +
    geom_line(color = "purple") +
    geom_point(color = "purple", size = 1) +
    labs(
      title = "Lasso Performance",
      x = "Number of Features (genes)",
      y = "Accuracy (Cross-Fold Validation)"
    ) +
    theme_pretty()
  
  lasso_perf
  
  ggsave(sprintf("Figures/lasso_perf_%s.png", s), 
         plot = lasso_perf,
         width = 2000,
         height = 1000,
         units = "px")
  
  lasso_tuning <- ggplot(lasso_cv_data, aes(x = -1*log(lambda), y = cvm)) +
    geom_point(color = "purple", size = 1) +
    geom_line(color = "purple", size = 0.8) +
    geom_ribbon(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), 
                fill = "purple", alpha = 0.2) +
    geom_vline(xintercept = -1*log(lasso.fit$lambda.min), 
               linetype = "dashed", color = "black") +
    geom_vline(xintercept = -1*log(lasso.fit$lambda.1se), 
               linetype = "dashed", color = "black") +
    scale_x_continuous(
      sec.axis = sec_axis(trans = ~., 
                          breaks = -1*log(lasso.fit$lambda[seq(1, length(lasso.fit$lambda), by = 9)]),
                          labels = colSums(lasso.fit$glmnet.fit$beta != 0)[seq(1, length(lasso.fit$lambda), by = 9)],
                          name = "Number of Features")
    ) +
    labs(
      title = "LASSO Tuning Curve",
      x = "-Log(λ)",
      y = "Binomial Deviance"
    ) + theme_pretty()
  
  lasso_tuning
  
  ggsave(sprintf("Figures/lasso_tuning_%s.png", s), 
         plot = lasso_tuning,
         width = 2000,
         height = 1000,
         units = "px")
  
  # END OF PLOTS
  
  
  # Extract nonzero coefficients (important genes)
  coeffs <- coef(lasso.fit, s = "lambda.1se")
  .coeffs_df <- as.data.frame(as.matrix(coeffs))
  
  # create lasso_genes df
  lasso_genes <- .coeffs_df |>
    filter(lambda.1se != 0)
  lasso_genes <- lasso_genes |>
    filter(rownames(lasso_genes) != "(Intercept)") |>
    arrange(desc(abs(lambda.1se)))
  
  # write to CSV, in order
  # rownames are the gene names
  write.csv(lasso_genes, file=sprintf("Data/key_genes_ranked_lasso_%s.csv", s))
  
  # Write names to text file
  cat(rownames(lasso_genes), file = sprintf("Data/key_genes_lasso_%s.txt", s))
  
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
  
  # PRETTYPLOT - Histogram lasso.pred ---
  lassopred_hist <- ggplot(lasso.pred, aes(x = lambda.1se)) + 
    geom_histogram(fill = "purple") +
    labs(
      title = "Distribution of LASSO coefficients",
      x = "Coefficient value (lambda.1se)"
    ) +
    theme_pretty()
  
  lassopred_hist
  
  ggsave(sprintf("Figures/lassopred_hist_%s.png", s), 
         plot = lassopred_hist,
         width = 2000,
         height = 1000,
         units = "px")
  # END OF PLOT ---
  
  cm_lasso <- confusionMatrix(pred_class, y_test)
  
  # Check ROC
  roc_obj <- roc(y_test, as.numeric(lasso.pred))
  auc(roc_obj)
  plot(roc_obj)
  
  # PRETTYPLOT - Lasso ROC ---
  roc_data <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  roc_lasso <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_path(size = 1, color = "purple") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#898989") +
    labs(
      title = "Lasso ROC Curve",
      x = "False Positive Rate",
      y = "True Positive Rate",
      subtitle = paste("AUC =", round(auc(roc_obj), 3))
    ) +
    coord_equal() +
    theme_pretty()
  
  roc_lasso
  
  ggsave(sprintf("Figures/roc_lasso_%s.png", s), 
         plot = roc_lasso,
         width = 1000,
         height = 1000,
         units = "px")
  # END OF PLOT ---
  
  
  ###### WARNING - GENAI! ----------------------------------
  # Used for permutation testing
  
  if (.permute){
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
  }
}

#-----------------------------------------------------------------------------
# Random forest
#-----------------------------------------------------------------------------
if (.RF){
  print("Running RF")
  control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
  
  RF.model <- train(x = x_train,
                    y = y_train,
                    method = "ranger", # for variable importance metrics
                    trControl = control,
                    ntree = 100,
                    importance = "permutation")
  
  RF.model
  
  # PRETTYPLOT - RF Performance ---
  RF_perf <- RF.model$results |> 
    ggplot(aes(x = mtry, y = Accuracy, color = splitrule, group = splitrule)) +
    geom_line() +
    geom_point() +
    labs(title = "RF Performance",
         x = "Number of Variable Splits",
         y = "Accuracy (Cross-Fold Validation)") +
    theme_pretty()
  
  RF_perf
  
  ggsave(sprintf("Figures/RF_perf_%s.png", s), 
         plot = RF_perf,
         width = 2000,
         height = 1000,
         units = "px")
  # END OF PLOT ---
  
  # Evaluate performance on test set
  # To get probability scores
  RF.prob <- predict(RF.model, newdata = x_test, type = "prob")
  
  # Threshold 0.5 showed best performance
  pred_class <- ifelse(RF.prob$NF > 0.5, "NF", "DCM")
  pred_class <- factor(pred_class, levels = levels(y_test))
  
  cm_RF <- confusionMatrix(pred_class, y_test)
  
  hist(RF.prob$DCM)
  hist(RF.prob$NF)
  
  # Check ROC plot
  roc_obj <- roc(y_test, as.numeric(pred_class))
  auc(roc_obj)
  plot(roc_obj)
  
  # PRETTYPLOT - ROC RF ---
  roc_data <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  roc_RF <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_path(size = 0.5, color = "purple") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    labs(
      title = "RF ROC Curve",
      x = "False Positive Rate",
      y = "True Positive Rate",
      subtitle = paste("AUC =", round(auc(roc_obj), 3))
    ) +
    coord_equal() +
    theme_pretty()
  
  roc_RF
  
  ggsave(sprintf("Figures/roc_RF_%s.png", s), 
         plot = roc_RF,
         width = 1000,
         height = 1000,
         units = "px")
  # END OF PLOT ---
  
  # Export genes considered most important 
  # Picked 38 to export since LASSO has 38-gene list
  
  key_genes_all <- varImp(RF.model)$importance |> arrange(desc(Overall))
  View(key_genes_all) 
  write.csv(key_genes_all,sprintf("Data/genes_RF_ranked_%s.csv", s))
  
  gene_names <- rownames(key_genes_all)
  cat(gene_names[1:38], file = sprintf("Data/key_genes_RF_%s.txt", s))
  cat(gene_names, file = sprintf("Data/key_genes_ALL_RF_%s.txt", s))
  
  ###### WARNING - GENAI! ----------------------------------
  if (.permute){
    # Try with randomly shuffled labels
    y_train_rand <- sample(y_train)
    RF.model_rand <- train(x = x_train,
                           y = y_train,
                           method = "ranger", # for variable importance metrics
                           trControl = control,
                           ntree = 100,
                           importance = "permutation")
    RF.model_rand$results
  }
}


#-----------------------------------------------------------------------------
# SVM-RFE
# Source used for help:
# https://www.geeksforgeeks.org/machine-learning/svm-feature-selection-in-r-with-example/
#-----------------------------------------------------------------------------
if (.SVM){
  print("Running SVM-RFE")
  control <- rfeControl(functions = caretFuncs, 
                        method = "cv", 
                        number = 10,
                        verbose = TRUE)
  
  svm_rfe <- rfe(x = x_train, 
                 y = y_train,
                 sizes = c(1:38),  # Number of features to select
                 preProc = c("center","scale"),
                 rfeControl = control,
                 method = "svmLinear")
  
  # Filter out the version using all genes
  svm_vars <- svm_rfe$variables |> 
    filter(Variables != 20781)
  
  # Get gene list in order of importance
  svm_genes_ranked <- svm_vars |> 
    group_by(var) |> 
    summarise(
      mean_Overall = mean(Overall),
      frequency = n()
    ) |> 
    arrange(desc(mean_Overall))
  
  svm_genes_ranked
  
  write.csv(svm_genes_ranked, sprintf("Data/svm-rfe_genes_ranked_%s.csv", s))
  
  # PRETTYPLOT - SVM-RFE Performance ---
  svm_perf <- svm_rfe$results |> 
    filter(Variables != 20781) |> 
    ggplot(aes(x = Variables, y = Accuracy)) +
    geom_line(color = "purple") +
    geom_point(color = "purple") +
    labs(title = "SVM-RFE Performance",
         x = "Number of Features (genes)",
         y = "Accuracy (Cross-Fold Validation)") +
    theme_pretty()
  
  svm_perf
  
  ggsave(sprintf("Figures/svmrfe_perf_%s.png", s), 
         plot = svm_perf,
         width = 2000,
         height = 1000,
         units = "px")
  # END OF PLOT ---
  
  SVMRFE.pred <- predict(svm_rfe, x_test)
  
  cm_SVMRFE <- confusionMatrix(SVMRFE.pred, y_test)
  
  # Check ROC plot
  roc_obj <- roc(y_test, as.numeric(SVMRFE.pred))
  auc(roc_obj)
  plot(roc_obj)
  
  # PRETTYPLOT - SVM RFE ROC ---
  roc_data <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  roc_SVM <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_path(size = 1, color = "purple") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    labs(
      title = "SVM-RFE ROC curve",
      x = "False Positive Rate",
      y = "True Positive Rate",
      subtitle = paste("AUC =", round(auc(roc_obj), 3))
    ) +
    coord_equal() +
    theme_pretty()
  
  roc_SVM
  
  ggsave(sprintf("Figures/roc_SVM-RFE_%s.png", s), 
         plot = roc_SVM,
         width = 1000,
         height = 1000,
         units = "px")
  
  # END OF PLOT ---
}



#-----------------------------------------------------------------------------
# Comparing performance
# 
#   Plots etc
#-----------------------------------------------------------------------------

# COMBINED PLT ---
if (.combine){
  combined_ROC <- roc_lasso + roc_SVM 
  ggsave(sprintf("Figures/combined_ROCs_%s.png", s), 
         plot = combined_ROC,
         width = 2000,
         height = 1000,
         units = "px")
  # END ---
  
  lasso_table <- as.data.frame(cm_lasso$overall)
  SVMRFE_table <- as.data.frame(cm_SVMRFE$overall)
  
  
  summary_stats <- data.frame(
    lasso = lasso_table,
    SVM_RFE = SVMRFE_table
  )
  
  summary_stats <- as.data.frame(t(summary_stats))
  
  write.csv(summary_stats, sprintf("Data/ML_summarystats_%s.csv", s))
}

