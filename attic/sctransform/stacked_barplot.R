library(Seurat)
library(Tweedieverse)
library(tidyverse)
library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(Maaslin2)
library(Tweedieverse)
library(omicsArt)
library(scanpy)

library(Seurat)
library(patchwork)
library(ggplot2)
install.packages("scCustomize")
library("scCustomize")
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )+
    ylab(feature) +
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 45, hjust=1)
          plot.margin = plot.margin, )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



features <- c("Col1a1", "Fbln2", "Adgre1", "Ptprc", "Itga8", "Tnnt3", "Pcdh9", "Flt1",
              "Ebf1", "Nrros", "Auts2", "Tshz2", "Myh11","Mecom","Spag17"
              ,"Ttn","Slc1a3","Cmss1","Gpc5","Far2","St6galnac3")

# Subset data.frame
#p1<-StackedVlnPlot(seur_obj,features) 
#p1+ theme(axis.text.y = element_text(size = 0))  + theme(axis.text.x = element_text(angle = 45, hjust=1))




p1 <- Stacked_VlnPlot(seur_obj,features,x_lab_rotate = TRUE) 

p1 +  theme(axis.text.y = element_text(size = 0))
