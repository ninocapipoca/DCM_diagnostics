#=============================================================================#
# projectPlots.R                                                              #
#	                                                                            #
# Creation of any additional plots we might want to include in the manuscript #
# that require information from all (or lots of) corners of our analysis      #
#                                                                             #
# NOTE - uses output CSV files from DataPreparationExploratoryAnalysis.R      #
#=============================================================================#

.packages <- c("ggplot2", "dplyr")
lapply(.packages, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

# Import data, filtered and imputed
gxData <- read.csv("Data/gxData_female_corr.csv", as.is = T, row.names = 1)
metadata <- read.csv("Data/metadata_female_corr.csv", as.is = T, row.names = 1)

key_genes <- scan("Data/AAA_COMMON_GENES.txt", what = "character")

#-----------------------------------------------------------------------------#
# Define theme
#-----------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------#
# Key genes plot(s)
#-----------------------------------------------------------------------------#

kg_plotdata <- as.data.frame(t(gxData[key_genes, ])) |>
  cbind(etiology = metadata[rownames(kg_plotdata), "etiology"]) |>
  pivot_longer(cols = -etiology, names_to = "Gene", values_to = "Expression")

kg_plot <- ggplot(kg_plotdata, aes(x = etiology, y = Expression, fill = etiology)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_pretty() +
  labs(title = "Gene Expression by Disease Status", 
       x = "",
       y = "Expression")

kg_plot

ggsave("Figures/keygene_boxplots.png", plot = kg_plot, width = 12, height = 6, dpi = 300)
