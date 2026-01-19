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
.compare <- TRUE # load both datasets for comparison

if (.compare){
  netw_genes_m <- read.csv("Data/geneInfo_males_top50.csv")
  lasso_genes_m <- read.csv("Data/key_genes_ranked_lasso_M.csv")
  SVM_genes_m <- read.csv("Data/svm-rfe_genes_ranked_M.csv")
  
  netw_genes_f <- read.csv("Data/geneInfo_top50.csv")
  lasso_genes_f <- read.csv("Data/key_genes_ranked_lasso_F.csv")
  SVM_genes_f <- read.csv("Data/svm-rfe_genes_ranked_F.csv")
  
  # Some formatting
  names(lasso_genes_m)[names(lasso_genes_m)=="X"] <- "Gene.ID"
  names(lasso_genes_f)[names(lasso_genes_f)=="X"] <- "Gene.ID"
  
  names(SVM_genes_m)[names(SVM_genes_m)=="var"] <- "Gene.ID"
  names(SVM_genes_f)[names(SVM_genes_f)=="var"] <- "Gene.ID"
} else {
  lasso_genes <- read.csv(sprintf("Data/key_genes_ranked_lasso_%s.csv", s))
  SVM_genes <- read.csv(sprintf("Data/svm-rfe_genes_ranked_%s.csv", s))
  
  if (s == "F"){
    netw_genes <- read.csv("Data/geneInfo_top50.csv")
  } else {
    netw_genes <- read.csv("Data/geneInfo_males_top50.csv")
  }
}


#-----------------------------------------------------------------------------#
# Venn Diagram (individual)
#-----------------------------------------------------------------------------#
if (!.compare){
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
    filename = sprintf("Figures/VennDiagram_genes_%s.png", s),
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
  cat(common_genes, file = sprintf("Data/AAA_COMMON_GENES_%s.txt", s))
}

#-----------------------------------------------------------------------------#
# Comparing M and F overlaps
#-----------------------------------------------------------------------------#

if (.compare) {
  mlist = list(SVM_genes_m, netw_genes_m, lasso_genes_m)
  flist = list(SVM_genes_f, netw_genes_f, lasso_genes_f)
  namelist = c("SVM-RFE", "WGCNA", "LASSO")
  
  for (i in 1:3){
    .male_data <- mlist[[i]]
    .female_data <- flist[[i]]
    print(i)
    cutoff <- min(nrow(.male_data), nrow(.female_data))
    
    genes <- list(male = .male_data$Gene.ID,
                  female = .female_data$Gene.ID)
    
    #genes <- lapply(genes, function(x) x[1:cutoff])
    
    venn.diagram(
      x = genes,
      category.names = c("Male", "Female"),
      filename = sprintf("Figures/VennDiagram_%s_comp.png", namelist[i]),
      output = TRUE,
      
      # Text
      # Numbers
      cex = 1.5,
      fontfamily = "sans",
      
      # Colors
      fill = c("#E8B4D4", "#B4E8D4"),  # Soft pink and soft teal
      alpha = c(0.5, 0.5),
      lty = 'blank',
      
      # Titles
      main = namelist[i],
      main.cex = 3,
      main.fontface = "bold",
      
      cat.cex = 1.8,
      cat.fontface = "bold",
      cat.pos = c(-27, 27),
      cat.dist = c(0.055, 0.055),
      cat.fontfamily = "sans",
    )
  }
  
  
  
}
