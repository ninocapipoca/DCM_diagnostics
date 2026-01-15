# -------------------------------------------------------------
#                                                             |                                                                 
# preprocess_ML.R                                             |
# Very simple data preparation to create ML-ready data        |
# NOTE - no imputation yet to avoid data leakage              |
#                                                             |
# -------------------------------------------------------------

# Load packages
library(dplyr)

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
# Filter dataset
#-----------------------------------------------------------------------------

# Log2 normalize
data <- log2(data + 1)

# Filter out SVA columns
metadata <- metadata_SVA |> dplyr::select(1:16)


.female <- TRUE # Flag for testing, if necessary
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

# Save as csv
write.csv(gxData, "Data/ML_gxData.csv")
write.csv(metadata, "Data/ML_metadata.csv")
