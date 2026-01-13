#=============================================================================#
# DataPreparationExploratoryAnalysis.R                                        #
#	                                                                            #
# Descriptive analysis, generation of participant summary table, PCA          #
#                                                                             #
# Start date: 08/01/2026      									                              #
#=============================================================================#

# If you are reading this, you are in the version of this file that uses
# the batch-corrected dataset from the MAGNet GitHub repo
# Link: https://github.com/mpmorley/MAGNet/tree/master

# Load packages ------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(gt) 
library(gtsummary)
library(pcaMethods)
library(limma)
library(sva)

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

# # Uncomment this code to redownload dataset
# # Requires CPMS_SVA_corrected.RDS from the GitHub above
# write.csv(CPMS_SVA_corrected,"Data/CPMS_SVA_corrected.csv")

# Load gene expression data and participant information (metadata)
# NOTE - Gene expression dataset contains log2-transformed CPM
gxData_all <- read.csv("Data/CPMS_SVA_corrected.csv", as.is = T, row.names = 1)
metadata_all <- read.csv("MAGNET_GX_2025/MAGNET_SampleData_18112022.csv", as.is = T, row.names = 1)

# Filter out only female participants who either have DCM or are healthy
metadata <- metadata_all |> 
  filter(gender == "Female" & (etiology == "DCM" | etiology == "NF"))

# Filter gxData_all to include only participants also in metadata
gxData <- gxData_all[, colnames(gxData_all) %in% rownames(metadata)]

# Turn all character-type columns into factors (categorical variables)
metadata <- metadata |> mutate_if(is.character, as.factor)

# Export filtered gxData and metadata for use in other scripts
write.csv(gxData, "Data/gxData_female_corr.csv")
write.csv(metadata, "Data/metadata_female_corr.csv")

# Participant table -----------------------------------------------------------
participant_table <- metadata |>
  # Define table structure and variables 
  tbl_summary(
    missing = "no",
    include = c('race', 'Diabetes', 
                'afib', 'Hypertension'),
    type = all_continuous() ~ "continuous2",
    by = etiology,
    label = list(afib = "Afibrillation", 
                 race = "Race"),
    # Def stats to be included; non-missing val, median, quartiles, range
    statistic = list(all_continuous() ~ c(
      "{N_nonmiss}",
      "{median} ({p25}, {p75})",
      "{min}, {max}"
    ),
    # Def stats for categorical; count & percent total
    all_categorical() ~ "{n} ({p}%)"
    )) |>
  # Statistical comparison accross groups ("is diff btwn them significant?")
  add_p() |>
  # Improve readability and information in table
  modify_header(label ~ "**Variable**") |>
  #modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Condition**") |>
  #modify_caption("**Table 1. Patient Characteristics**<br>
  #               Patient characteristics grouped by type of cardiomyopathy: dilated (DCM),
  #               hypertrophic (HCM), peripartum (PPCM) and healthy individuals (NF)") |>
  bold_labels() |>
  # Save as png
  as_gt() |>
  gtsave("participant_table.png", path=paste0(directory_path, '/Figures'))

# For display in Rstudio
participant_table


#-----------------------------------------------------------------------------#
# Diagnostic plots
#-----------------------------------------------------------------------------#

# Box & density plots --------------------------------------------------------

# Sanity check to make sure dimensions are compatible
all(rownames(metadata) == colnames(gxData)) 

# Create plots for each condition
for (cond in c("NF", "DCM")) {
  plot_data <- gather(gxData[, metadata$etiology == cond], 
                      key = "SampleID", value = "CPM")
  print(cond)
  print(dim(plot_data))
  
  # Box plot for each condition
  # if (length(unique(plot_data$SampleID)) > 200) {
  #   # If too many samples, becomes unreadable; plot first 60
  #   box_plot <- plot_data |>
  #     filter(SampleID %in% unique(SampleID)[1:60]) |>
  #     ggplot(aes(x = SampleID, y = CPM)) + 
  #     geom_boxplot() +
  #     theme_minimal(base_size = 7) +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #     theme(axis.text.x = element_blank()) +
  #     labs(title = paste("Gene expression accross random samples -", cond),
  #          x = "Sample (label removed)",
  #          y = "Log2 transformed CPM")
  # }
  
#  else {
    box_plot <- plot_data |>
      ggplot(aes(x = SampleID, y = CPM)) + 
      geom_boxplot() +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Gene expression accross random samples -", cond),
           x = "Sample ID",
           y = "Log2 transformed CPM")
#  }
  
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
  
  ggsave(sprintf("Density plot for %s Samples.png", cond), plot = density_plot, 
         path=paste0(directory_path, '/Figures'))
  
  ggsave(sprintf("Gene Expression for %s Samples.png", cond), plot = box_plot, 
         path=paste0(directory_path, '/Figures'))
  
}

# PCA Plot(s) ----------------------------------------------------------------

# nipalsPCA chosen instead of default svd bcs of missing vals - note: slower.
pca_res <- pca(t(gxData), nPcs = 10, method = "nipals")

# Sanity check - see if sample names correspond
all(rownames(pca_res@scores) == rownames(metadata))

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
ggsave("Scree_plot.png", plot = scree, path=paste0(directory_path, '/Figures'))

# Compare PCs; see if clustering arises corresponding to any metadata feature
plot_pca <- cbind(data.frame(pca_res@scores), metadata)


pc1_label = paste("PC1 (score ", round(var_expl[1,2], digits=2), ")", sep="")
pc2_label = paste("PC2 (score ", round(var_expl[2,2], digits=2), ")", sep="")

# Scatterplots
# Split suggests PC2 may be related to etiology
pca1v2 <- ggplot(plot_pca, aes(x = PC1, y = PC2)) + 
  geom_point(aes(col = etiology)) +
  labs(x = pc1_label, y = pc2_label)
  theme_minimal()

pca1v2
ggsave("PC1_PC2_corr.png", plot = pca1v2, path=paste0(directory_path, '/Figures'))

#-----------------------------------------------------------------------------#
# Differential Expression Analysis
#-----------------------------------------------------------------------------#

design <- model.matrix(~ 0 + etiology+age+race, data = metadata)
fit <- lmFit(gxData, design)

# Contrast DCM and control
cont.matrix <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge_res_DCM <- topTable(ebFit, coef = 'DCMvsControl', number = nrow(gxData))

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

DCM_diffexp <- get_sig(dge_res_DCM)

#-----------------------------------------------------------------------------#
# Pathway Enrichment Analysis
#-----------------------------------------------------------------------------#

# Get up & downregulated genes
# Up & downregulated genes for each contrast
DCM_upreg <- DCM_diffexp[DCM_diffexp$logFC > 0, ]
DCM_downreg <- DCM_diffexp[DCM_diffexp$logFC < 0, ]

# As a precaution, create copies
univ <- DCM_diffexp
up <- DCM_upreg
down <- DCM_downreg

# Add gene_id columns to each dataframe
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

# KEGG Enrichment analyis ----------------------------------------------------

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
write.csv(kegg_down@result, "Data/kegg_down_DCM_corr.csv")
write.csv(kegg_up@result,"Data/kegg_up_DCM_corr.csv")

# Create dot plots showing most enriched categories and how much
dotplot_up <- dotplot(kegg_up, showCategory = 10) +
  ggtitle("KEGG PEA - Upregulated") +
  theme(axis.text.y = element_text(size = 8))
ggsave("Dotplot_Upreg_KEGG_corr.png", path=paste0(directory_path, '/Figures'))

dotplot_down <- dotplot(kegg_down, showCategory = 10) +
  ggtitle("KEGG PEA - Downregulated") +
  theme(axis.text.y = element_text(size = 8))
ggsave("Dotplot_Downreg_KEGG_corr.png", path=paste0(directory_path, '/Figures'))

# GO Enrichment Analysis -----------------------------------------------------

GO_up <- enrichGO(gene = up_new$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'ENTREZID',
                      ont = 'All',
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

GO_down <- enrichGO(gene = down_new$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENTREZID',
                    ont = 'All',
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

# Save output as CSV
write.csv(GO_up@result, "Data/GO_up_DCM_corr.csv")
write.csv(GO_down@result,"Data/GO_down_DCM_corr.csv")

# Create GO dot plots
dotplot_up <- dotplot(GO_up, showCategory = 10) +
  ggtitle("GO PEA - Upregulated") +
  theme(axis.text.y = element_text(size = 8))
ggsave("Dotplot_Upreg_GO_corr.png", path=paste0(directory_path, '/Figures'))

dotplot_down <- dotplot(kegg_down, showCategory = 10) +
  ggtitle("GO PEA - Downregulated") +
  theme(axis.text.y = element_text(size = 8))
ggsave("Dotplot_Downreg_GO_corr.png", path=paste0(directory_path, '/Figures'))
