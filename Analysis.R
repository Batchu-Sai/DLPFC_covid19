
# Script to create figures of SARS-COV2 relevant genes taken from Figure 1 (https://www.nature.com/articles/s41586-021-03710-0#Fig1)
# Data taken from https://doi.org/10.1038/s41593-020-00787-0

rm(list = ls())
library(cowplot)
library(ggplot2)
library(stringr)
library(spatialLIBD)
library(RColorBrewer)
library(viridis)
source('Utils.R') # Code modified for figures, original code from spatialLIBD
setwd("~/Research/COVID19_Spatial_DLPFC/Analysis")
options(stringsAsFactors=FALSE)

# Download data into folder

spe <- fetch_data(type = 'spe', destdir = "~/Research/COVID19_Spatial_DLPFC/Data")

# Choosing the last four samples, have all layers as seen in original manuscript
# First plot the layers of these samples in one column

samples <- spe@int_metadata[["imgData"]]@listData[["sample_id"]]

layer_plots <- list()

for (i in 1:length(samples[(9:12)])){
  image1 <- vis_clus2(
    spe = spe,
    clustervar = "spatialLIBD",
    alpha = 0.45,
    point_size = 0.85,
    sampleid = as.character(samples[[i]]),
    colors = libd_layer_colors)+
    ggplot2::labs(title = paste0("Sample ", as.character(samples[[i]])))+
    theme(plot.title = element_text(size=14, face="bold"))
  layer_plots[[i]] <- image1
  
}

layer_plot <- plot_grid(plotlist=layer_plots, ncol = 1)
png("~/Research/COVID19_Spatial_DLPFC/Figures/Layer_Plot.png", width = 450, height = 1050)
print(layer_plot)
dev.off()



### Plot genes of interest

# View genes available
d <- as.matrix(rowData(spe)$gene_search)

genes_interest <- c("ACE2; ENSG00000130234",
                    "BSG; ENSG00000172270", 
                    "NRP1; ENSG00000099250", 
                    "NRP2; ENSG00000118257",
                    
                    "TMPRSS2; ENSG00000184012",
                    "TMPRSS11A; ENSG00000187054",
                    "FURIN; ENSG00000140564",
                    "CTSB; ENSG00000164733",
                    "CTSL; ENSG00000135047",
                    
                    "LY6E; ENSG00000160932",
                    "IFITM1; ENSG00000185885",
                    "IFITM2; ENSG00000185201",
                    "IFITM3; ENSG00000142089",
                    "IFNAR1; ENSG00000142166",
                    "IFNAR2; ENSG00000159110"
                    )

gene_plots <- list()  

for (g in genes_interest){
  for (i in 1:length(samples[(9:12)])){
    image1 <- vis_gene2(
      spe,
      sampleid = as.character(samples[(9:12)][[i]]),
      geneid = g,
      spatial = TRUE,
      assayname = "logcounts",
      minCount = 0,
      viridis = F,
      image_id = "lowres",
      cont_colors = c("aquamarine4","springgreen", "goldenrod", "red"),
      alpha = 0.45,
      point_size = 1.2,
      ... = ""
    )+theme(plot.title = element_text(size=16, face="bold"))
    gene_plots[[i]] <- image1
    plot_final <- plot_grid(plotlist=gene_plots, ncol = 1)
    png(paste0("~/Research/COVID19_Spatial_DLPFC/Figures/", g, ".png"), width = 450, height = 1050)
    print(plot_final)
    dev.off()
  }
}
