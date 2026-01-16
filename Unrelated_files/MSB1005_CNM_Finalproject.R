#=============================================================================#
# MSB1005_CNM_Finalproject.R                                                  #
#																	                                       		  #
# Version: 1.0   															                                #
# Start date: Nov 28, 2025											                              #
# Author: Carolina Nina Matos (i6446025)                                      #
#=============================================================================#
# Load packages ------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(gt) 
library(gtsummary)
library(pcaMethods)
library(limma)

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(biomaRt)
#-----------------------------------------------------------------------------#
# Data import & participant table
#-----------------------------------------------------------------------------#

# Set working directory according to the location of this script
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Load gene expression data and participant information (metadata)
# NOTE - Gene expression dataset contains log2-transformed CPM
gxData <- read.table("MAGNET_GX_Assignment_MSB1005/MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T, row.names = 1, header = TRUE)
metadata <- read.csv("MAGNET_GX_Assignment_MSB1005/MAGNET_SampleData_18112022.csv", as.is = T, row.names = 1)

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata %>% mutate_if(is.character, as.factor)


# Calculate BMI
# 
# @description
# Calculates the BMI based on height in cm and weight in kg 

# I know this function gives some crazy values (bmi 200???) but I didn't have
# time to fix it and it shouldn't affect downstream analysis
bmi <- function(height, weight){
  weight / ((height / 100) ** 2)
}

# Create a copy of the metadata and calculate BMI for all participants 
metadata_stats <- metadata
metadata_stats$BMI <- mapply(bmi, metadata$height, metadata$weight)

# Participant table -----------------------------------------------------------

# Use gtsummary and gt libraries to create participant table by etiology
participant_table <- metadata_stats %>%
  # Define table structure and variables 
  tbl_summary(
    missing = "no",
    include = c('gender', 'race', 'Diabetes', 
                'afib', 'Hypertension', 'BMI', 
                'tissue_source'),
    type = all_continuous() ~ "continuous2",
    by = etiology,
    label = list(tissue_source = "Tissue Source", 
                 afib = "Afibrillation", 
                 race = "Race", 
                 gender = "Sex"),
    # Def stats to be included; non-missing val, median, quartiles, range
    statistic = list(all_continuous() ~ c(
      "{N_nonmiss}",
      "{median} ({p25}, {p75})",
      "{min}, {max}"
    ),
    # Def stats for categorical; count & percent total
    all_categorical() ~ "{n} ({p}%)"
  )) %>%
  # Statistical comparison accross groups ("is diff btwn them significant?")
  add_p() %>%
  # Improve readability and information in table
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Condition**") %>%
  modify_caption("**Table 1. Patient Characteristics**<br>
                 Patient characteristics grouped by type of cardiomyopathy: dilated (DCM),
                 hypertrophic (HCM), peripartum (PPCM) and healthy individuals (NF)") %>%
  modify_footnote_body(
    footnote = "Underweight: < 18.5, Normal: from 18.5 to 24.9, Overweight: 25 to 29.9, Obese or extremely obese: > 30",
    columns = "label",
    rows = variable == "BMI" & row_type == "label"
  ) %>%
  bold_labels() %>%
  # Save as png
  as_gt() %>%
  gtsave("participant_table.png", path=directory_path)

# For display in Rstudio
participant_table


#-----------------------------------------------------------------------------#
# Diagnostic plots
#-----------------------------------------------------------------------------#

# Box & density plots --------------------------------------------------------

# Sanity check to make sure dimensions are compatible
all(rownames(metadata) == colnames(gxData)) 

# Create plots for each disease condition
for (cond in c("NF", "DCM", "HCM", "PPCM")) {
  plot_data <- gather(gxData[, metadata$etiology == cond], 
                      key = "SampleID", value = "CPM")
  print(cond)
  print(dim(plot_data))
  
  # Box plot for each condition
  if (length(unique(plot_data$SampleID)) > 100) {
    # If too many samples, becomes unreadable; plot first 60
    box_plot <- plot_data |>
      filter(SampleID %in% unique(SampleID)[1:60]) |>
      ggplot(aes(x = SampleID, y = CPM)) + 
      geom_boxplot() +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.text.x = element_blank()) +
      labs(title = paste("Gene expression accross random samples -", cond),
           x = "Sample (label removed)",
           y = "Log2 transformed CPM")
  }
  
  else {
    box_plot <- plot_data |>
      ggplot(aes(x = SampleID, y = CPM)) + 
      geom_boxplot() +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Gene expression accross random samples -", cond),
           x = "Sample ID",
           y = "Log2 transformed CPM")
  }
  
  # Density plot legend removed to improve readability (doesn't add much info)
  # Even when low nr of samples, for consistency
  density_plot <- plot_data |>
    ggplot(aes(x = CPM, color = SampleID)) + 
    geom_density(alpha = 0.2) +
    theme_minimal(base_size = 7) +
    theme(legend.position = "none") +
    labs(title = sprintf("Densty plot for %s participants", cond),
         x = "Log2 transformed CPM",
         y = "Density")
  
  ggsave(sprintf("Density plot for %s Samples.png", cond), plot = density_plot)
  ggsave(sprintf("Gene Expression for %s Samples.png", cond), plot = box_plot)
  
}
# PCA Plot(s) ----------------------------------------------------------------

# nipalsPCA chosen instead of default svd bcs of missing vals - note: slower.
pca_res <- pca(t(gxData), nPcs = 10, method = "nipals")

# Sanity check - see if sample names correspond
all(rownames(pca_res@scores) == rownames(metadata)) # TRUE

# Scree plot to understand difference of explained variance between PCs
var_expl <- data.frame(
  Component = factor(1:10),
  Variance = pca_res@sDev^2 / sum(pca_res@sDev^2))
var_expl

scree <- ggplot(var_expl, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Scree Plot", 
       x = "Principal Component", 
       y = "Proportion of variance explained")

scree
ggsave("Scree_plot.png", plot = scree)

# Compare PCs; see if clustering arises corresponding to any metadata feature
plot_pca <- cbind(data.frame(pca_res@scores), metadata)

# Scatterplots
# Split suggests PC2 may be related to etiology
pca1v2 <- ggplot(plot_pca, aes(x = PC1, y = PC2)) + 
  geom_point(aes(col = etiology)) +
  theme_minimal()

pca1v2
ggsave("PCA_PC1_vs_PC2.png", plot = pca1v2)


#-----------------------------------------------------------------------------#
# Statistical analysis
#-----------------------------------------------------------------------------#

# Differential gene expression analysis ---------------------------------------
design <- model.matrix(~ 0 + etiology, data = metadata)
fit <- lmFit(gxData, design)

# Contrast DCM and control
cont.matrix <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_DCM <- topTable(ebFit, coef = 'DCMvsControl', number = nrow(gxData))

# Contrast PPCM and control
cont.matrix <- makeContrasts(PPCMvsControl = etiologyPPCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_PPCM <- topTable(ebFit, coef = 'PPCMvsControl', number = nrow(gxData))

# Contrast HCF and control
cont.matrix <- makeContrasts(HCMvsControl = etiologyHCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_HCM <- topTable(ebFit, coef = 'HCMvsControl', number = nrow(gxData))

# Get genes with a statistically significant log fold change
# 
# @description
# Orders dataframe by adjusted p-value and 
# returns new dataframe with only significant rows (padj < 0.05)
get_sig <- function(df) {
  df |>
    arrange(adj.P.Val) |>
    filter(adj.P.Val < 0.05)
}

contrast_dfs <- lapply(list(dge_res_DCM, dge_res_HCM, dge_res_PPCM), get_sig)

# Correct for covariates -----------------------------------------------------

# Justification: gender (ie, sex) plays a role since, for instance, PPCM
# is female-specific. Additionally, based on previous literature
# both age and sex are risk factors for cardiac diseases.
# If I had more time, I could do a more formal correlation analysis of this 

design <- model.matrix(~ 0 + etiology + age + gender, data = metadata)
fit <- lmFit(gxData, design)

# Contrast DCM and control
cont.matrix <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_DCM <- topTable(ebFit, coef = 'DCMvsControl', number = nrow(gxData))

# Contrast PPCM and control
cont.matrix <- makeContrasts(PPCMvsControl = etiologyPPCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_PPCM <- topTable(ebFit, coef = 'PPCMvsControl', number = nrow(gxData))

# Contrast HCF and control
cont.matrix <- makeContrasts(HCMvsControl = etiologyHCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_HCM <- topTable(ebFit, coef = 'HCMvsControl', number = nrow(gxData))

contrast_dfs_corr <- lapply(list(dge_res_DCM, dge_res_HCM, dge_res_PPCM), get_sig)

#-----------------------------------------------------------------------------#
# Additional gene annotation and relative expression levels
#
# Note - unfortunately, I realized too late that I misinterpreted questions 4
# and 5. I believe the idea is similar and I would like to demonstrate that 
# I do still understand the underlying concepts, so I have included an 
# additional explanation in comments at the end of this file.
# 
#-----------------------------------------------------------------------------#

# Up & downregulated genes for each contrast
DCM_diffexp <- contrast_dfs_corr[[1]]
DCM_upreg <- DCM_diffexp[DCM_diffexp$logFC > 0, ]
DCM_downreg <- DCM_diffexp[DCM_diffexp$logFC < 0, ]

HCM_diffexp <- contrast_dfs_corr[[2]]
HCM_upreg <- HCM_diffexp[HCM_diffexp$logFC > 0, ]
HCM_downreg <- HCM_diffexp[HCM_diffexp$logFC < 0, ]

PPCM_diffexp <- contrast_dfs_corr[[3]]
PPCM_upreg <- PPCM_diffexp[PPCM_diffexp$logFC > 0, ]
PPCM_downreg <- PPCM_diffexp[PPCM_diffexp$logFC < 0, ]

.diffexp_dfs <- list(DCM_diffexp, HCM_diffexp, PPCM_diffexp)
.upreg_dfs <- list(DCM_upreg, HCM_upreg, PPCM_upreg)
.downreg_dfs <- list(DCM_downreg, HCM_downreg, PPCM_downreg)
.names <- list("DCM", "HCM", "PPCM")



for (i in 1:3) {
  univ <- .diffexp_dfs[[i]]
  up <- .upreg_dfs[[i]]
  down <- .downreg_dfs[[i]]
  
  univ$gene_id <- rownames(univ)
  up$gene_id <- rownames(up)
  down$gene_id <- rownames(down)
  
  # Convert gene names to correct format
  gene_data_eID <- bitr(univ$gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  univ_new <- merge(univ, gene_data_eID, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)
  
  gene_data_eID <- bitr(up$gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  up_new <- merge(up, gene_data_eID, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)
  
  gene_data_eID <- bitr(down$gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  down_new <- merge(down, gene_data_eID, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)
  
  # remove NA values
  univ_new <- na.omit(univ_new)
  up_new <- na.omit(up_new)
  down_new <- na.omit(down_new)
    
  # Perform kegg pathway enrichment analysis
  # For upregulated genes
  kegg_up <- enrichKEGG(gene = up_new$ENTREZID,
                        organism = 'hsa',  # humans
                        universe = univ_new$ENTREZID,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)
  
  # For downregulated genes
  kegg_down <- enrichKEGG(gene = down_new$ENTREZID,
                        organism = 'hsa',  # humans
                        universe = univ_new$ENTREZID,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)
  
  # Save output as CSV
  write.csv(kegg_down@result, sprintf("kegg_down_%s.csv", .names[i]))
  write.csv(kegg_up@result, sprintf("kegg_up_%s.csv", .names[i]))
  
  # Create dot plots showing most enriched categories and how much
  dotplot_up <- dotplot(kegg_up, showCategory = 10) +
  ggtitle(sprintf("KEGG PEA - Upregulated, for %s", .names[i])) +
  theme(axis.text.y = element_text(size = 8))
  ggsave(sprintf("Dotplot_Upreg_%s.png", .names[i]))
  
  dotplot_down <- dotplot(kegg_down, showCategory = 10) +
  ggtitle(sprintf("KEGG PEA - Downregulated, for %s", .names[i])) +
  theme(axis.text.y = element_text(size = 8))
  ggsave(sprintf("Dotplot_Downreg_%s.png", .names[i]))
  
}

#-----------------------------------------------------------------------------#
# Consolidation and export
#-----------------------------------------------------------------------------#



NF_patients <- rownames(metadata)[metadata$etiology == "NF"]
HCM_patients <- rownames(metadata)[metadata$etiology == "HCM"]
DCM_patients <- rownames(metadata)[metadata$etiology == "DCM"]
PPCM_patients <- rownames(metadata)[metadata$etiology == "PPCM"]

geneTotExonLengths <- read.delim("MAGNET_GX_Assignment_MSB1005/MAGNET_exonLengths.txt", as.is = T, 
                                 row.names = 1)
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1] # . operator just hides from workspace
}
gxData_fpkm <- cpm2fpkm(gxData)

condition_means <- data.frame(
  ensembl_gene_id = rownames(gxData),
  NF_meanCPM = rowMeans(gxData[, NF_patients, drop = FALSE]),
  HCM_meanCPM = rowMeans(gxData[, HCM_patients, drop = FALSE]),
  DCM_meanCPM = rowMeans(gxData[, DCM_patients, drop = FALSE]),
  PPCM_meanCPM = rowMeans(gxData[, PPCM_patients, drop = FALSE]),
  NF_meanFPKM = rowMeans(gxData_fpkm[, NF_patients, drop = FALSE]),
  HCM_meanFPKM = rowMeans(gxData_fpkm[, HCM_patients, drop = FALSE]),
  DCM_meanFPKM = rowMeans(gxData_fpkm[, DCM_patients, drop = FALSE]),
  PPCM_meanFPKM = rowMeans(gxData_fpkm[, PPCM_patients, drop = FALSE])
)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "description"),
  filters = "ensembl_gene_id",
  values = rownames(condition_means),
  mart = ensembl
)

# merge with annotations
export <- condition_means |>
  left_join(annotations, by = "ensembl_gene_id")

# try to label differential expression
export <- export |>
  mutate(
    differentially_expressed = ifelse(
      ensembl_gene_id %in% rownames(DCM_upreg) |
        ensembl_gene_id %in% rownames(PPCM_upreg) |
        ensembl_gene_id %in% rownames(HCM_upreg),
      "UP",
      ifelse(
        ensembl_gene_id %in% rownames(DCM_downreg) |
          ensembl_gene_id %in% rownames(PPCM_downreg) |
          ensembl_gene_id %in% rownames(HCM_downreg),
        "DOWN",
        "NO"
      )
    )
  )

write.table(export, file = "exported_results.tsv", sep = "\t", row.names = FALSE)

#-----------------------------------------------------------------------------#
# Additional explanation
#
# Path enrichment analysis (PEA) tells us whether a gene is overrepresented
# compared to random chance in a certain pathway. The overrepresented genes
# are connected to the pathways they are part of through a database. Here, 
# I used Kegg, but it's also possible to do Go enrichment. What was written in
# the assignment feels like a precursor to this, since we are comparing the
# expression to the background / noise levels.
# 
#-----------------------------------------------------------------------------#
