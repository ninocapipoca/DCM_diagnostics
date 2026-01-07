#=============================================================================#
# MSB1005 Final Assignment                                                    #
#																	                                       		  #
#                                                                             #
# Date: 12 Dec 2025																	                          #														                                              #   
#=============================================================================#




#=============================================================================#
#1) Data import:                                                              #
#=============================================================================#

# a) Importing all the data files
setwd(".") 
# Importing Gene expression dataset which is a table which contains all the genes
gxData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T,
                     row.names = 1)
# Importing sample information file which contains the clinical data of the patients
sampleInfo <- read.csv("MAGNET_SampleData_18112022.csv", as.is = T, 
                       row.names = 1)
# Importing the csv file which contains the exon length for every gene ensemble ID
geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T,
                         row.names = 1)


# b) Export a publication-ready table of participant characteristics, including 
#statistics comparing the four etiologies.
require(gtsummary)
# Select the following fields for the publication-ready table:
#   etiology - used to classify participant characteristics into four categories
#   race - relevant for patients genetics 
#   gender - important biological variable, especially for the PPCM disease (female only)
#   lv_mass - measure of heart enlargement, important when measuring cardiomyopathy
#   VTVF - ventricular tachycardia/ventricular fibrillation,showing arrhythmia complications
#   Diabetes - relevant condition when considering overall characteristics
#   Hypertension - major cardiovascular risk factor
Publication_table <- sampleInfo %>% 
  select(etiology, race, age, gender, lv_mass, VTVF, Diabetes, Hypertension) %>%
  tbl_summary(by = etiology, missing = "no") %>% # Detects variable types and applies appropriate statistics
  add_overall() %>%  #Adds extra column showing statistics for all etiologies combined 
  add_p() # Performs appropriate statistical tests, comparing categorical continuous variables
write.csv(as_tibble(Publication_table), "summary.csv", row.names = FALSE) 


#=============================================================================#
#2) Diagnostic plots:                                                         #
#=============================================================================#

# a)At least one data distribution figure (e.g. boxplots, density plots) that enables 
#comparing samples.
require(tidyr)
require(ggplot2)
# Four boxplots are created (one per etiology group), showing the distribution of gene expression data within each etiology
# comparing to the healthy patients data
#First, a subset is created to extract NF columns and then gather is used to reshape the result to a different 
#format (sample IDs as factor). Each box represents the expression distribution 
#across all genes for that individual sample.
plotData <- gather(gxData[, sampleInfo$etiology == "NF"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()

# DCM patient samples. 
# Same subsetting and reshaping process.
plotData <- gather(gxData[, sampleInfo$etiology == "DCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()

# PPCM patient samples. 
# Same subsetting and reshaping process.
plotData <- gather(gxData[, sampleInfo$etiology == "PPCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()

# HCM patient samples. 
# Same subsetting and reshaping process.
plotData <- gather(gxData[, sampleInfo$etiology == "HCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()



# b) At least one PCA figure showing the sample clustering colored by relevant co-variates
#comparing samples.
require(pcaMethods)
pcaRes <- pca(t(gxData), nPcs = 10)
plot(pcaRes)
plotPcs(pcaRes, c(1,2))

# To create one relevant PCA figure, the following factors and relationships were chosen:
all(rownames(pcaRes@scores) == sampleInfo[, 1]) 
plotData <- cbind(data.frame(pcaRes@scores), sampleInfo)
# Plot two PCs to capture variance. Point size scaled by age and color by disease etiology. 
# Age - etiology checks if either age or etiology contribute to the primary variation in gene expression
ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = age, col = etiology))
# Plot two PCs to capture variance. Point size scaled by age and color by race 
# Age - race plots the population spread with regards to race and age 
ggplot(plotData, aes(x = PC2, y = PC3)) + 
  geom_point(aes(size = age, col = race))
# Plot two PCs to capture variance. Point size scaled by lv_mass and color by Hypertension 
# lv_mass - Hypertension is useful for checking biological variation when considering disease effects
ggplot(plotData, aes(x = PC4, y = PC5)) + 
  geom_point(aes(size = lv_mass, col = Hypertension))
# Plot two PCs to capture variance. Point size scaled by LVEF and color by Diabetes. 
# LVEF - Diabetes checks if the two components reveal a pattern when assessing coexistent conditions 
ggplot(plotData, aes(x = PC6, y = PC7)) + 
  geom_point(aes(size = LVEF, col = Diabetes))
# Plot two PCs to capture variance. Point size scaled by LVEF and color by Diabetes
# LVEF - Hypertension (same as the previous relationship, but with hypertension)
ggplot(plotData, aes(x = PC6, y = PC8)) + 
  geom_point(aes(size = LVEF, col = Hypertension))
# weight - Diabetes checks if these two components are responsible for the remaining variation not explained by the main components
# Plot two PCs to capture variance. Point size scaled by weight and color by Diabetes 
ggplot(plotData, aes(x = PC9, y = PC10)) + 
  geom_point(aes(size = weight, col = Diabetes))


#=============================================================================#
#3) Statistical analysis:                                                     #
#=============================================================================#

# a) Perform a differential gene expression analysis comparing DCM, HCM and 
# PPCM patients to the healthy donors.
# b) Correct for relevant co-variates and add comments to the scripts explaining 
#your choices. 
require(limma)
# The following covariates were chosen:
# - age: accounts for age-related gene expression changes
# - gender: responsible for gender-specific expression differences
# - race: adjusts for population stratification effects
# - TIN.median: corrects for RNA sample quality
# Build a matrix with etiology as the main factor including covariates (age, gender, race, TIN).
# Fits linear model to expression data and defines the contrast 
# between the two main factors contrast and applies
# moderation to stabilize variance estimates and extracts the results.
design_DCM <- model.matrix(~0 +etiology + age + gender + race + TIN.median., data = sampleInfo)
fit_DCM <- lmFit(gxData, design_DCM)
cont.matrix_DCM <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF,
                             levels = design_DCM)
fit_DCM <- contrasts.fit(fit_DCM, cont.matrix_DCM)
ebFit_DCM <- eBayes(fit_DCM, trend =TRUE)
dgeRes_DCM <- topTable(ebFit_DCM, coef = 'DCMvsControl', number = nrow(gxData))

# The same process as the previous etiology
# Since PPCM can only take place in females (postpartum), gender is excluded in this case
design_PPCM <- model.matrix(~0 +etiology + age + race + TIN.median., data = sampleInfo)
fit_PPCM <- lmFit(gxData, design_PPCM)
cont.matrix_PPCM <- makeContrasts(PPCMvsControl = etiologyPPCM - etiologyNF,
                             levels = design_PPCM)
fit_PPCM <- contrasts.fit(fit_PPCM, cont.matrix_PPCM)
ebFit_PPCM <- eBayes(fit_PPCM, trend =TRUE)
dgeRes_PPCM <- topTable(ebFit_PPCM, coef = 'PPCMvsControl', number = nrow(gxData))

# The same process as the previous two etiologies
# HCM affects both genders, thus gender is included
design_HCM <- model.matrix(~0 +etiology + age + gender + race + TIN.median., data = sampleInfo)
fit_HCM <- lmFit(gxData, design_HCM)
cont.matrix_HCM <- makeContrasts(HCMvsControl = etiologyHCM - etiologyNF,
                             levels = design_HCM)
fit_HCM <- contrasts.fit(fit_HCM, cont.matrix_HCM)
ebFit_HCM <- eBayes(fit_HCM, trend =TRUE)
dgeRes_HCM <- topTable(ebFit_HCM, coef = 'HCMvsControl', number = nrow(gxData))


#=============================================================================#
# 4) Additional gene annotation:                                              #
#=============================================================================#

# a) Retrieve gene symbols and gene names based on the provided Ensembl gene 
#identifiers. 
# Retrieve gene notations, symbols, description and chromosome names (needed for the noise calculation later) for all genes in gxData using biomaRt 
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(gxData)
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id','hgnc_symbol', 'description', 'chromosome_name'),
  filters = 'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

# b) Merge this additional annotation with the gene expression data object.
# Merge gene annotations with expression data, converting gxData rownames
# to a column, then join with gene_annotations table on matching IDs
gxData_annotated <- merge(
  data.frame(ensembl_id = rownames(gxData), gxData),
  gene_annotations,
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE  
)


#=============================================================================#
#5) Relative expression levels:                                               #
#=============================================================================#

# a)Transform the data to FPKM values
# Use the function to convert both gxData and gxData_annotated to FPKM values
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1]
}
gxData_fpkm <- cpm2fpkm(gxData)



#b) Assess for each gene in the dataset whether it is expressed above 
# background (noise) level. 
# For example: you can base this on comparing the average expression of a 
# gene to the average expression of Y chromosome genes in female subjects. 

# Check which columns contain numeric data
numeric_cols <- sapply(gxData_annotated, is.numeric)
string_cols <- !numeric_cols

# For FPKM conversion, the numeric columns are extracted from gxData_annotated
gxData_numeric <- gxData_annotated[, numeric_cols]
gxData_string <- gxData_annotated[, string_cols, drop = FALSE]

# Only the numeric columns are converted to FPKM
gxData_fpkm <- cpm2fpkm(gxData_numeric)

# Obtain the IDs of all female participants from sample info 
female_patients <- rownames(sampleInfo[sampleInfo$gender == "Female", ])

# Keep only the IDs which appear in gxData_fpkm
female_patients <- female_patients[female_patients %in% colnames(gxData_fpkm)]

# Extract the IDs of the genes on the Y chromosome
y_chromosome_genes <- rownames(gxData_string)[gxData_string$chromosome_name == "Y"]

# Again, only keep the ones which appear in gxData_fpkm
y_chromosome_genes <- y_chromosome_genes[y_chromosome_genes %in% rownames(gxData_fpkm)]

# Calculate the noise level of Y chromosome genes in females
noise_expression <- mean(as.matrix(gxData_fpkm[y_chromosome_genes, female_patients]), na.rm = TRUE)

# Calculate the mean FPKM data for each gene
gene.means_fpkm <- rowMeans(gxData_fpkm, na.rm = TRUE)

# Combine the separated numeric and string data back into one structure
gxData_annotated_fpkm <- cbind(gxData_string, gxData_fpkm)
rownames(gxData_annotated_fpkm) <- gxData_annotated_fpkm$ensembl_id
gxData_annotated_fpkm$mean_fpkm <- gene.means_fpkm

#Establish which genes are expressed above noise level
gxData_annotated_fpkm$expressed <- gene.means_fpkm > noise_expression



#=============================================================================#
# 6) Export the results:                                                      #
#=============================================================================#

# Extract the patient IDs from sampleInfo and only keep 
# those which exist in gxData for all etiologies 

DCM_samples <- rownames(sampleInfo[sampleInfo$etiology == "DCM", ])
DCM_samples <- DCM_samples[DCM_samples %in% colnames(gxData)]

HCM_samples <- rownames(sampleInfo[sampleInfo$etiology == "HCM", ])
HCM_samples <- HCM_samples[HCM_samples %in% colnames(gxData)]

PPCM_samples <- rownames(sampleInfo[sampleInfo$etiology == "PPCM", ])
PPCM_samples <- PPCM_samples[PPCM_samples %in% colnames(gxData)]

NF_samples <- rownames(sampleInfo[sampleInfo$etiology == "NF", ])
NF_samples <- NF_samples[NF_samples %in% colnames(gxData)]
all_genes <- rownames(gxData)
results <- data.frame(ensembl_id = all_genes, row.names = all_genes)

# Add annotation data and match by rownames
annotation_cols <- gxData_annotated_fpkm[, c(names(gxData_string), "mean_fpkm", "expressed")]
results <- merge(results, annotation_cols, 
                          by = "row.names", all.x = TRUE)
rownames(results) <- results$Row.names
results$Row.names <- NULL
results$ensembl_id <- NULL 

# Calculate mean expression for each group
results$mean_expr_DCM <- rowMeans(gxData[rownames(results), DCM_samples], na.rm = TRUE)
results$mean_expr_HCM <- rowMeans(gxData[rownames(results), HCM_samples], na.rm = TRUE)
results$mean_expr_PPCM <- rowMeans(gxData[rownames(results), PPCM_samples], na.rm = TRUE)
results$mean_expr_Control <- rowMeans(gxData[rownames(results), NF_samples], na.rm = TRUE)

# Add DGEA results for all etiologies
dgeRes_DCM_subset <- dgeRes_DCM[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
colnames(dgeRes_DCM_subset) <- paste0("DCM_", colnames(dgeRes_DCM_subset))
results <- merge(results, dgeRes_DCM_subset, 
                          by = "row.names", all.x = TRUE)
rownames(results) <- results$Row.names
results$Row.names <- NULL

dgeRes_HCM_subset <- dgeRes_HCM[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
colnames(dgeRes_HCM_subset) <- paste0("HCM_", colnames(dgeRes_HCM_subset))
results <- merge(results, dgeRes_HCM_subset, 
                          by = "row.names", all.x = TRUE)
rownames(results) <- results$Row.names
results$Row.names <- NULL

dgeRes_PPCM_subset <- dgeRes_PPCM[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
colnames(dgeRes_PPCM_subset) <- paste0("PPCM_", colnames(dgeRes_PPCM_subset))
results <- merge(results, dgeRes_PPCM_subset, 
                          by = "row.names", all.x = TRUE)
rownames(results) <- results$Row.names
results$Row.names <- NULL

results$ensembl_id.y <- NULL  
results$ensembl_id.x <- NULL


# Export all the combined results to one file
write.table(results, 
            file = "Final_results.txt", 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA, 
            quote = FALSE)


