.packages <- c("dplyr", "pROC", "caret")
lapply(.packages, require, character.only = TRUE)


# Set working directory
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load data
tryCatch({
  # Try loading ML-ready gxData
  data <- read.csv("Data/ML_gxData.csv", row.names = 1, header = TRUE)
}, warning = function(w) {
  cat("A warning occurred:", conditionMessage(w), "\n")
  print("If the file ML_gxData.csv does not exist, run preprocess_ML.R first")
})

tryCatch({
  # Try loading ML-ready metadata
  metadata_SVA <- read.csv("Data/ML_metadata.csv", row.names = 1, header = TRUE)
}, warning = function(w) {
  cat("A warning occurred:", conditionMessage(w), "\n")
  print("Try downloading the file & adding it to the Data folder: https://www.dropbox.com/s/eihem5fbnkg7bpm/phenoData.csv?dl=0")
})
