#=================================================#
# MSB1005 R Assignment                            #
#																	                #
# Alma C Verme, i6265381		                      #
# Date: Dec 12, 2025                              #
#=================================================#

# Disclaimer: I am confused about how much you expect us to explain and have gone for the safe bet: a 
# lot (sorry! if not necessary). Since it is an assessment, I felt I had to explain myself more than
# I would have in a typical portable script. I quickly started running out of time at the end, and accepted rushing.

# First we load all the necessary packages for this script:
library(tidyverse)
library(gtsummary)
library(broom)
library(webshot2)
library(pcaMethods)
library(limma)
library(biomaRt)

#==================================================================================================#
# 1a. Data import
#     Set your working directory and import the two data files. To import gxData, we use read.table() 
#     with arguments 'header = TRUE' to tell R that the first row in the data is a header. To import 
#     sampleInfo we use read.csv and use the 'row.names = 1' for the first column to be the row titles.

setwd("/Users/mikiverme/Downloads/MAGNET_GX_2025") # <-- Set your own working directory here!

gxData <- read.table("MAGNET_GeneExpressionData_CPM_19112020.txt", header = TRUE)
sampleInfo <- read.csv("MAGNET_SampleData_18112022.csv") 
#     Reminder: the "sample_name" column from the .csv file is now the row names. It is no longer
#     longer a variable + observations in the data.

#-----------------------------------------------------------------------------#
# 1b. Export publication-ready table
#     I realised we'd need a package for this task and googled: "packages to make publication-ready 
#     tables in R" to get an idea, for which I found the overwhelming majority (in Reddit/other forums) 
#     to use gtsummary. 
#     See: https://r-graph-gallery.com/367-automate-data-summary-table-with-gtsummary.html

#     To automate the creation of the table, start by making a copy of the sampleInfo data (to always 
#     keep one clean version). 

sampleInfo.with_BMI <- sampleInfo

#     Already knowing the data, we go straight to selecting which variables to include. I chose to 
#     show the patient BMI (calculated from height and weight), despite both variables individually 
#     having significant differences across the four etiology groups  (p-value = 0.046 and 0.002 
#     respectively). To avoid redundancy the table only displays the BMI. The dplyr package was used 
#     to mutate the BMI column to the new dataframe.

sampleInfo.with_BMI <- sampleInfo.with_BMI %>% 
  dplyr::mutate(BMI = weight/(height/100)^2)
print(sampleInfo.with_BMI$BMI)   # quick check to see if BMI values make sense and are in the dataframe

#     Additionally, we do not include some variables, eg., minexpr since it is the minimum expression 
#     of a gene across all patients and will therefore show no variation when tested. To create the 
#     table we use the tbl_summary() function, with the arguments 'by' determining which variables to 
#     classify the information by and 'include = c()' determining which variables to include. The 
#     'label' argument is used to title the variables appropriately for publication (those that need 
#     it). The table creation automatisation includes the structure by default (ie., 'n(%)' and 
#     'Median (Q1,Q3)').

sampleInfo.publish <- sampleInfo.with_BMI %>%
  tbl_summary(by = etiology, include = c(tissue_source, gender, age, race, BMI, hw, lv_mass, afib, VTVF, Diabetes, Hypertension, RIN, TIN.median.),
              label = list(tissue_source = "Tissue source",
                        gender = "Gender",
                        age = "Age (years)",
                        race = "Race",
                        hw = "Heart weight (g)",
                        lv_mass = "Left ventricle mass (g)",
                        afib = "Atrial fibrilation",
                        VTVF = "Ventricular tachycardia/fibrilation",
                        TIN.median. = "Median TIN value"
                        )) %>%
  bold_labels() %>%  # makes the labels bold to distinguish variables from their (sub)levels
  add_overall() %>%
  add_p()  # Runs one of the significance tests per variable, depending on variable type

#     If the variable is categorical, the function runs Fisher's exact test which does not assume 
#     normality and for continuous variables the Kruskal-Wallis rank sum test is performed, also 
#     doesn't assume normality (is non-parametric). A normality test could be run before making the 
#     table to observe the behaviour of the data and then a more specific test could be defined in 
#     this add_p() argument. On the table itself, the 'Unknown' rows  identify 
#     the number of 'NA's found in that variable. Though it adds a lot more visual information to 
#     the table, I think it is good to keep it for transparency. To remove the unknown, set the 
#     following argument for the tbl_summary() function: 'missing = "no"'.

#     To export the table, we use gtsummary's function gtsave() to save it in the current working 
#     directory folder. I chose to save it as a PNG file, but this can be substituted by .jpg (JPG), 
#     .htm (HTML), etc. depending on the preference.

sampleInfo.publish %>% 
  as_gt() %>%   # converts the gt_summary data into a gt table object
  gt::gtsave("sampleInfo.png")   # saves the gt table as a PNG file 

#==================================================================================================#
# 2a. Diagnostic plots
#     Diagnostic plots let you visualise the distribution of the data. First, we must rearrange the
#     data using the pivot_longer() function so that every gene expression is a sample (row). For every 
#     row, we then have the gene ID, a patient ID, and their respective expression of that gene. 
#     Additionally, to be able to separate the data by etiology, we use the left_join() function. 
#     Because it expects a dataframe, but when I did it the first time I thought it was not useful to 
#     add ALL the sampleInfo columns to the long-format gxData. Therefore within the left_join function 
#     we use select() to choose only sample_name and etiology, and merge them with gxData by patientID.

all(rownames(sampleInfo) == colnames(gxData))

gxData2 <- gxData %>%
  pivot_longer(cols = -EnsemblGeneID,   # Makes an individual row for every gene ID per patients and respective expression value
    names_to = "SampleID",   # previously the column names, patients are now listed under a column called SampleID
    values_to = "CPM")   # expression values aligned with respective gene and patient, column is called CPM
gxData2 <- gxData2 %>%
  left_join(sampleInfo %>% select(sample_name, etiology), by = c("SampleID" = "sample_name"))  

#     Then we can visualise the data with a boxplot, having each patient along the x-axis and their 
#     gene expression (CPM) in the y-axis. I chose to put each etiology as a separate panel on the 
#     same figure (named: 'boxplot'), but plots can be done separately if you're looking for better 
#     visibility. Though, the patient ID labels will never be visible with so many samples. The last 
#     graph is a combined graph of the four individual ones if one would want it (recommended to be 
#     viewed large). A density plot can also be useful to observe the distribution of all expression.

boxplot <- ggplot(gxData2, aes(SampleID, CPM, color = etiology)) +
  geom_boxplot() +
  facet_wrap(~ etiology)   # creates a separate figure panel for each etiology

density <- ggplot(gxData2, aes(x = CPM, color = etiology)) +
  geom_density()

#     If you want to have a closer look at the box plot information, you can create individual plots separated
#     by etiology. Below is an example, subsetting for NF, that you could copy-paste and fill in according to
#     the rest of the etiologies.

NF.plot <- gxData2 %>%   # Title according to etiology chosen
  filter(etiology == "NF") %>%   # choose which etiology to filter TRUE by
  ggplot(aes(x = SampleID, y = CPM)) +
  geom_boxplot(color = "turquoise") +  # choose the colour to distinguish one graph from another
  labs(title = "NF")   # titles the plot

#-----------------------------------------------------------------------------#
# 2b. PCA
#     NUMBER principal component analyses are performed to visualise potential covariates. These include
#     every sampleInfo column that is not the sample name itself and their respective etiology. First, we
#     must create a new gxData.num since we need PCA data to be numeric and we previously chose to keep 
#     the EnsembleGeneID column (when importing the dataset). 
#     Additionally, the pcaMethods package expects samples to be rows, so we use the pca() function on the 
#     gxData transposed, and set it to perform the analysis on the first 10 components.

gxData.num <- gxData %>% 
  column_to_rownames("EnsemblGeneID") # make a dataframe where the EnsemblGeneID are the row names

pca <- pca(t(gxData.num), nPcs = 10) # set pca() to analyse the transposed data and set number of components to analyse
pca@R2  # check the PCA scores (make sure they're decreasing, PC1 should have the highest weight/value)

#     To visualise confounding effects on the arrangement of samples across the PCA components, we have to
#     save the PCA analysis data (scores) and merge them with the sampleInfo data to plot them together.

pca.scores <- as.data.frame(pca@scores)  # save PCA scores as a dataframe
pca.plotData <- cbind(pca.scores, sampleInfo[, -1])  # bind pca.scores with sampleInfo (ignoring its first patient ID column, 
                                                     # since it's pca.scores' rownames)

#     Then we can plot the components against each other and colour by gender (or other confounders):

pca.plot1 <- ggplot(pca.plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(col = etiology, size = age))

pca.plot2 <- ggplot(pca.plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(col = gender, shape = etiology))

#-----------------------------------------------------------------------------#
# 2c. Export plots
#     The two final plots are exported as PNGs, the combined boxplot and density plot. Below them, there
#     is an example (reminder, rather!) for saving the boxplots separated by etiology, were you to choose 
#     to use them. 

ggsave("combined_boxplot.png", plot = boxplot)  # Combined boxplot
ggsave("density.png", plot = density)  # Density plot
ggsave("pca_plot1.png", plot = pca.plot1)  # PCA plot 1 (etiology, age) # there is a clear difference between the controls and the patients
ggsave("pca_plot2.png", plot = pca.plot2)  # PCA plot 2 (etiology, gender)

# Separate etiology boxplot reminder (2a):
ggsave("NF_plot.png",  # <-- Title the plot and identify the saving format after the period
       plot = NF.plot)   # <-- Choose which plot you want to save by filling the argument

#==================================================================================================#
# 3a. Differential Gene Expression Analysis
#     gxData must be numeric, so we will be using the gxData.num dataframe for all of these questions. We
#     also have to make etiology a factor, to classify the data by DCM, HCM, PPCM or NF (no other possible
#     answers). 

sampleInfo$etiology <- factor(sampleInfo$etiology, levels = c("NF","DCM","HCM","PPCM"))

#     Now we make a design matrix, with 0 separating the etiologies (there's no intercept between them),
#     fit a linear regression for each gene, and also define the 'contrasts' (the comparisons that the 
#     model aims to analyse).

design.matrix <- model.matrix(~ 0 + etiology, data = sampleInfo)
fit <- lmFit(gxData.num, design.matrix)  # this fits a linear regression for each gene
cont.matrix <- makeContrasts(
  DCMvsNF  = etiologyDCM - etiologyNF, 
  HCMvsNF  = etiologyHCM - etiologyNF,
  PPCMvsNF = etiologyPPCM - etiologyNF,
  levels = design.matrix)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge.DCM <- topTable(ebFit, coef = 'DCMvsNF', number = nrow(gxData.num))
dge.HCM <- topTable(ebFit, coef = 'HCMvsNF', number = nrow(gxData.num))
dge.PPCM <- topTable(ebFit, coef = 'PPCMvsNF', number = nrow(gxData.num))

#     I used the head() function to inspect the results, which are a dataframe containing the log fold 
#     change (gene up/down-regulation compared to the control), its avg. expression, the t-stat (cal-
#     culated earlier by the eBayes() function), the p-value, a p-value adjusted for multiple testing 
#     (false discovery rate), and  the odds that the gene is truly differentially expressed (log; the 
#     higher the value, the more likely it is) To inspect the significant results of the DGE we can use 
#     the following summaries:

summary(dge.DCM$adj.P.Val < 0.05)  # number of significant genes, DCM
summary(dge.HCM$adj.P.Val < 0.05)  # HCM
summary(dge.PPCM$adj.P.Val < 0.05)  # PPCM

#-----------------------------------------------------------------------------#
# 3b. Adjust for relevant covariates
#     We use the same name of the dataframe, so we don't have to run the contrasts again <-- don't forget
#     or get them confused!

design.matrix <- model.matrix(~ 0 + etiology + age + gender, data = sampleInfo) # Adding age and gender to adjust the model for
fit <- lmFit(gxData.num, design.matrix) 
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
dge.DCM <- topTable(ebFit, coef = 'DCMvsNF', number = nrow(gxData.num))
dge.HCM <- topTable(ebFit, coef = 'HCMvsNF', number = nrow(gxData.num))
dge.PPCM <- topTable(ebFit, coef = 'PPCMvsNF', number = nrow(gxData.num))

#==================================================================================================#
# 4a. Additional gene annotation
#     To retrieve gene symbols & names from Ensembl, we first have to retrieve from our data a list of the
#     gene IDs. Then, with the useEnsembl() function (connecting to the Ensembl database) we select "genes"
#     from the homosapiens Dataset. I spent way too long searching for what the attributes we wanted are,
#     but anyways we identify them as the gene ID and the HGNC symbol (nomeclature), from which we want to 
#     retrieve the specific ones we listed (with the rownames() function) and from the 'human.mart' mart.

genes_IDs <- rownames(gxData.num)

human.mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), # I added 'chromosome_name' because we will need it for 5b
                     values = genes_IDs,
                     mart = human.mart)

#-----------------------------------------------------------------------------#
# 4b. Merge Ensembl data
#     Because we left the original gxData with the first column as the EnsembleGeneID, we can now use that
#     dataframe directly. First, we change the first column's name to match that of the gene IDs in gxData. 
#     That way, we can merge the new annotations with gxData by gene ID, and with the argument 'all.y = TRUE' 
#     to keep all genes (rows) even if there is no annotation for them on Ensemble.

colnames(annotations)[1] <- "EnsemblGeneID" # Change the first column's name
gxData.ann <- merge(annotations, gxData, by = "EnsemblGeneID", all.y = TRUE) # merge annotations and gxData by gene ID

#==================================================================================================#
# 5a. Relative expression levels
#     In week 4 we created the function cpm2fpkm which converts the CPM values back to linear scale (from
#     log) and multiplies them by 1000 to get fpkm ("per kilobase" per patient; meaning we cannot use gene
#     data in this unit to compare samples). 
#     To start, we need the total exon lengths for each gene. *We use the numerical gxData.num again

geneTotExonLengths <- read.delim("4_MAGNET_exonLengths.txt", as.is = T, row.names = 1) # import from Week 4
all(rownames(geneTotExonLengths) == rownames(gxData.num)) # Check that they match!
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1]
}
gxData.fpkm <- cpm2fpkm(gxData.num) # use the conversion function

#-----------------------------------------------------------------------------#
# 5b. Relative expression levels
#     Analysing the female expression of Y-chromosome genes is apparently a great way to measure noise in
#     data (from own data; empirical and reducing a possible level of bias). First, one has to retrieve all 
#     patient IDs that are female and, separately, we retrieve all Y chromosome-coding genes.

females <- subset(sampleInfo, gender == "Female", select = sample_name) # subset females
genes.Y <- subset(annotations, chromosome_name == "Y", select = EnsemblGeneID) # subset Y-coding genes

#     Then we subset from the original gene expression dataset for Y-coding genes expressed by the females.
#     The 'drop = FALSE' argument forces the results to stay as a dataframe and not 'drop' to a vector. In 
#     the second line, we change the rownames to the EnsembleGeneID that matches the observation number.
#     Then, we can tell R to calculate a value for the average noise expression.

females.Y <- gxData[gxData[[1]] %in% genes.Y$EnsemblGeneID, colnames(gxData) %in% females$sample_name, drop = FALSE]
rownames(females.Y) <- gxData[gxData[[1]] %in% genes.Y$EnsemblGeneID, 1]
avg.gxNoise <- mean(as.matrix(females.Y))

#     We then have to calculate the average expression of every gene. By this point, I am very annoyed not to 
#     have imported gxData with rownames = 1 (I thought I was being smart). Anyways, we calculate the row means
#     in two steps (lines), and compare the average expression to the average noise. If the gene has more noise
#     than average it is 'noisy'.

genex <- colnames(gxData)[-1]
gxData$avg_expression <- rowMeans(gxData[, genex])
gxData$noisy <- gxData$avg_expression > avg.gxNoise

#==================================================================================================#
# 6a. Export the results
#     To export the results we merge the overall expression data that we have extracted for each gene
#     (average expression, whether it's noisy T/F, ) from gxData with gene symbol and chromosome name
#     from the annotations.

export.data <- merge(gxData[, c("EnsemblGeneID", "avg_expression", "noisy")],
                     annotations[, c("EnsemblGeneID", "hgnc_symbol", "chromosome_name")],
                     by = "EnsemblGeneID", all.x = TRUE)

#     We have to make sure that the CPM are numeric and calculate the mean, to compute the average 
#     expression per etiology.

gxData2$CPM <- as.numeric(as.character(gxData2$CPM)) # make sure that CPM is numeric to compute the mean
etiology.avg <- gxData2 %>%
  group_by(EnsemblGeneID, etiology) %>%  
  summarise(avg_expression_etiology = mean(CPM, na.rm = TRUE), .groups = "drop")

#     Then, we have to merge the two previous steps per gene.

export_data <- merge(export_data, etiology.avg, by = "EnsemblGeneID", all.x = TRUE)

#     Finally, to make the table more readible, we use pivot_wider() function with all the data we
#     we have extracted. We can use the write.table() function to export the data as a text table,
#     making sure it is delimited by tabs.

export_data_wide <- export_data %>%
  pivot_wider(id_cols = c(EnsemblGeneID, hgnc_symbol, chromosome_name, avg_expression, noisy),
    names_from = etiology, names_prefix = "avg_", values_from = avg_expression_etiology)  # use etiology names and prefix them with 'avg_' to remember!
write.table(export_data_wide, "GeneExpression_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(export_data_wide)
