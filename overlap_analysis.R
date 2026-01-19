#==============================================================================
# overlap_analysis.R                                                          #
#																	                                       		  #
# Analyzing the overlap between genes found by WGCNA, LASSO & SVM-RFE         #  															                                #
#=============================================================================#
.packages <- c("dplyr", "ggplot2", "VennDiagram", "RColorBrewer")
lapply(.packages, require, character.only = TRUE)


#-----------------------------------------------------------------------------#
# Data import
#-----------------------------------------------------------------------------#
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

s <- "F" # change this if needed

lasso_genes <- read.csv(sprintf("Data/key_genes_ranked_lasso_%s.csv", s))
SVM_genes <- read.csv(sprintf("Data/svm-rfe_genes_ranked_%s.csv", s))

if (s == "F"){
  netw_genes <- read.csv("Data/geneInfo_top50.csv")
} else {
  netw_genes <- read.csv("Data/geneInfo_males_top50.csv")
}


#-----------------------------------------------------------------------------#
# Venn Diagram
#-----------------------------------------------------------------------------#

# Make sure number of genes in each set is comparable
cutoff <- min(sapply(list(lasso_genes, SVM_genes, netw_genes), nrow))

genes <- list (lasso = lasso_genes$X,
               SVM = SVM_genes$var,
               netw = netw_genes$Gene.ID)

genes <- lapply(genes, function(x) x[1:cutoff])

venn.diagram(
  # Settings/structure here taken from
  # https://r-graph-gallery.com/14-venn-diagramm
  x = genes,
  category.names = c("SVM-RFE", "LASSO", "WGCNA"),
  filename = "Data/VennDiagram_genes.png",
  output = TRUE,
  
  # Colors!
  lwd = 2,
  lty = 'blank',
  fill = brewer.pal(3, "Pastel2"),
  
  # Text
  # Numbers
  cex = 1.5,
  fontfamily = "sans",
  
  # Titles
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Get common genes and check they're correct
common_genes <- Reduce(intersect, genes)

common_genes %in% genes$lasso
common_genes %in% genes$SVM
common_genes %in% genes$netw

# Hurray!! Export :)
cat(common_genes, file = "Data/AAA_COMMON_GENES.txt")