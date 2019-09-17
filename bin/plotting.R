library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)

date = Sys.Date()
source("bin/functions.R")

for (i in 1:4) {
  time_point <- list.files(file.path(getwd(), "data/timepoints"))[i]
  
  load(paste0("data/timepoints/",time_point))
  df <- data.frame('Cluster' = as.numeric(sc@cpart), sc@tsne) 
  df$ID <- names(sc@cpart)
  rownames(df) <- names(sc@cpart)
  df$Condition <- ifelse(grepl("HFD", df$ID), "HFD", "NC") %>%
    factor(levels = c("NC", "HFD"))
  
  #remove small clusters
  cell_numbers <- numeric()
  for (i in 1:max(df$Cluster, na.rm = T)) {
    cell_numbers[i] <- length(na.omit(df$Cluster[df$Cluster==i]))
  }
  names(cell_numbers) <- c(1:max(df$Cluster, na.rm = T))
  retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
  
  ord_clust <- clustheatmap(sc, final = T)
  
  ord_clust <- ord_clust[ord_clust %in% retain_cl]
  df$Cluster <- factor(df$Cluster, levels = ord_clust)
  
  df <- df[df$Cluster %in% retain_cl,]
  df$Cluster <- factor(df$Cluster, levels = ord_clust)
  colnames(df)[2:3] <- c("V1", "V2")
  
  #plot conditions
  tsne <- tsne_plot(df, FILL = df$Condition, fill_colors = colors_pat, point_outline = "black", point_size = 7, line_width = 0.25) +
    scale_fill_brewer(palette = "Set2")
  
  tsne
  
  ggsave(paste0('plots/tsne/', gsub("\\..*","", time_point), '-condition-tsne-plot.pdf'), width = 8.57, height = 5.79)  
  
  svg(paste0('plots/tsne/',  gsub("\\..*","", time_point), '-condition-tsne-plot.svg'), width = 8.57, height = 5.79)
  print(tsne)
  dev.off()
  
  #plot clusters
  tsne <- tsne_plot(df, FILL = df$Cluster, fill_colors = colors_many, point_outline = "black", point_size = 7, line_width = 0.25) 
  tsne
  
  ggsave(paste0('plots/tsne/', gsub("\\..*", "", time_point), '-clusters-tsne-plot.pdf'))  
  
  svg(paste0('plots/tsne/', gsub("\\..*","", time_point), '-clusters-tsne-plot.svg'), width = 8.57, height = 5.79)
  print(tsne)
  dev.off()
  
  print(i)
}


