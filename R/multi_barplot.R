library(tidyr)
library(dplyr)
library(reshape2)
devtools::install_github('omicsEye/deepath', force = TRUE)
library(deepath)
library("ggplot2")
library(cowplot)

install.packages(c('dplyr','pbapply', 'lme4', 'lmerTest', 'car', 'cplm', 'pscl', 'logging', 'ggrepel', 'gridExtra', 'future', 'cowplot', 'Hmisc', 'MultiDataSet', 'TSP', 'htmlTable', 'igraph', 'insight', 'lubridate', 'mgcv', 'mvtnorm', 'optparse', 'parameters', 'pillar', 'pkgload', 'plotly', 'rlang', 'rvest', 'seriation', 'usethis', 'viridis', 'xfun', 'yulab.utils', "labdsv", "seriation","diffusionMap"), repos='http://cran.r-project.org')
BiocManager::install("ropls")
devtools::install_github('omicsEye/omicsArt', force = TRUE)
install_github('omicsEye/omicsArt', force = TRUE)

library(omicsArt)
#setting the working directory
setwd("C:/Users/ranoj/Desktop/Single_Cell_output/analysis")

number_of_sig_to_keep <- 20
sig_threshold <- 0.05

## read for fibroblast
fibroblast <- read.delim(
  "matched/Tweedieverse_output_Fibroblasts/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
fibroblast_data_wound <- fibroblast[fibroblast$metadata=="Sample_type" & fibroblast$value=="Wound" ,]
rownames(fibroblast_data_wound) <- fibroblast_data_wound$feature


## read for macrophges
macrophages <- read.delim(
  "matched/Tweedieverse_output_Macrophages/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

macrophages_data_wound <- macrophages[macrophages$metadata=="Sample_type" & macrophages$value=="Wound" ,]
rownames(macrophages_data_wound) <- macrophages_data_wound$feature


## read for stromal
stromal <- read.delim(
  "matched/Tweedieverse_output_Stromal/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

stromal_data_wound <- stromal[stromal$metadata=="Sample_type" & stromal$value=="Wound" ,]
rownames(stromal_data_wound) <- stromal_data_wound$feature

##Using the defined genes which are markers for different cell types###
imp_genes <- c("CD68", "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1","Krt14", "Krt10", "Krt5","Plin1", "Adipoq", "Pparg", "Fabp4","Ptprc","Pecam1", "CD34")
#rownames(imp_genes) <- c("CD68", "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1","Krt14", "Krt10", "Krt5","Plin1", "Adipoq", "Pparg", "Fabp4","Ptprc","Pecam1", "CD34")

# use fibroblast as reference
order_sig <- rownames(fibroblast_data_wound)[1:number_of_sig_to_keep]
order_sig <- unique(append(order_sig,imp_genes))#Put a unique before append
##add the association for each celltype
fibroblast_data_wound <- fibroblast_data_wound[order_sig,]
##reinitalizing the rownames to the list 
rownames(fibroblast_data_wound) <- order_sig
fibroblast_data_wound[is.na(fibroblast_data_wound$coef),'coef'] = 0
fibroblast_data_wound[is.na(fibroblast_data_wound$pval),'pval'] = 1
fibroblast_data_wound[is.na(fibroblast_data_wound$qval),'qval'] = 1
fibroblast_data_wound$feature <- rownames(fibroblast_data_wound)
fibroblast_data_wound<- fibroblast_data_wound[order(fibroblast_data_wound$coef),]
order_sig_1 <- rownames(fibroblast_data_wound)
fibroblast_data_wound <- within(fibroblast_data_wound,
                                        feature <- factor(feature,
                                                          levels=order_sig_1))

macrophages_data_wound <- macrophages_data_wound[order_sig,]
rownames(macrophages_data_wound) <- order_sig
macrophages_data_wound[is.na(macrophages_data_wound$coef),'coef'] = 0
macrophages_data_wound[is.na(macrophages_data_wound$pval),'pval'] = 1
macrophages_data_wound[is.na(macrophages_data_wound$qval),'qval'] = 1
macrophages_data_wound$feature <- rownames(macrophages_data_wound)
macrophages_data_wound <- within(macrophages_data_wound,
                                            feature <- factor(feature,
                                                              levels=order_sig_1))


stromal_data_wound <- stromal_data_wound[order_sig,]
rownames(stromal_data_wound) <- order_sig
stromal_data_wound[is.na(stromal_data_wound$coef),'coef'] = 0
stromal_data_wound[is.na(stromal_data_wound$pval),'pval'] = 1
stromal_data_wound[is.na(stromal_data_wound$qval),'qval'] = 1
stromal_data_wound$feature <- rownames(stromal_data_wound)
stromal_data_wound <- within(stromal_data_wound,
                                           feature <- factor(feature,
                                                             levels=order_sig_1))

fibroblast_data_wound_temp_diff_bar <- diff_bar_plot(fibroblast_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
macrophages_data_wound_temp_diff_bar <- diff_bar_plot(macrophages_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')
stromal_data_wound_temp_diff_bar <- diff_bar_plot(stromal_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')



## read association
box_association <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-Omics/analysis/meatbolites_Tweedieverse/figures/Group_gg_associations.RDS")
## do plots

fig2_metabolites <- ggdraw() +
  draw_plot(fibroblast_data_wound_temp_diff_bar,
            x = 0, y = .47, width = .55, height = .53) +
  draw_plot(macrophages_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .55, y = .47, width = .225, height = .53) +
  draw_plot(stromal_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.line.y = element_blank()),
            x = .775, y = .47, width = .225, height = .53) +
  draw_plot(box_association[[4]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[15]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association_1[[14]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association_2[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a",  "Severe", "non-Severe", "non-COVID", "b", "c", "d", "e")),
                  size = 7,x = c(0, .28, .53, .76, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))
fig3_metabolites

ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)

box_association <- readRDS("matched/Tweedieverse_output_Fibroblasts/figures/Sample_type_gg_associations.RDS")
box_association_1 <- readRDS("matched/Tweedieverse_output_Macrophages/figures/Sample_type_gg_associations.RDS")
box_association_2 <- readRDS("matched/Tweedieverse_output_Stromal/figures/Sample_type_gg_associations.RDS")

fig_2 <- ggdraw() +
  draw_plot(fibroblast_data_wound_temp_diff_bar,
            x = 0, y = .47, width = .40, height = .55) +
  draw_plot(macrophages_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .40, y = .47, width = .35, height = .55) +
  draw_plot(stromal_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .75, y = .47, width = .25, height = .55)  +
  draw_plot(box_association[[4]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[15]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association_1[[14]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association_2[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a",  "Fibroblast", "Macrophages","Stromal" , "b", "c", "d", "e")),
                  size = 7,x = c(0, .10, .55, .80, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))

ggsave(filename = 'matched_barplot_fibroblast.pdf', plot=fig_2, width = 183, height = 110, units = "mm", dpi = 350)

### Epithelial as reference
fibroblast <- read.delim(
  "matched/Tweedieverse_output_Epithelial/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
fibroblast_data_wound <- fibroblast[fibroblast$metadata=="Sample_type" & fibroblast$value=="Wound" ,]
rownames(fibroblast_data_wound) <- fibroblast_data_wound$feature


## read for endothelial
macrophages <- read.delim(
  "matched/Tweedieverse_output_Endothelial/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

macrophages_data_wound <- macrophages[macrophages$metadata=="Sample_type" & macrophages$value=="Wound" ,]
rownames(macrophages_data_wound) <- macrophages_data_wound$feature


## read for Stem
stromal <- read.delim(
  "matched/Tweedieverse_output_Stem/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

stromal_data_wound <- stromal[stromal$metadata=="Sample_type" & stromal$value=="Wound" ,]
rownames(stromal_data_wound) <- stromal_data_wound$feature

imp_genes <- c("CD68", "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1","Krt14", "Krt10", "Krt5","Plin1", "Adipoq", "Pparg", "Fabp4","Ptprc","Pecam1", "CD34")

##Using the defined genes which are markers for different cell types###
order_sig <- rownames(fibroblast_data_wound)[1:number_of_sig_to_keep]
order_sig <- unique(append(order_sig,imp_genes))#Put a unique before append
##add the association for each celltype
fibroblast_data_wound <- fibroblast_data_wound[order_sig,]
##reinitalizing the rownames to the list 
rownames(fibroblast_data_wound) <- order_sig
fibroblast_data_wound[is.na(fibroblast_data_wound$coef),'coef'] = 0
fibroblast_data_wound[is.na(fibroblast_data_wound$pval),'pval'] = 1
fibroblast_data_wound[is.na(fibroblast_data_wound$qval),'qval'] = 1
fibroblast_data_wound$feature <- rownames(fibroblast_data_wound)
fibroblast_data_wound<- fibroblast_data_wound[order(fibroblast_data_wound$coef),]
order_sig_1 <- rownames(fibroblast_data_wound)
fibroblast_data_wound <- within(fibroblast_data_wound,
                                feature <- factor(feature,
                                                  levels=order_sig_1))

macrophages_data_wound <- macrophages_data_wound[order_sig,]
rownames(macrophages_data_wound) <- order_sig
macrophages_data_wound[is.na(macrophages_data_wound$coef),'coef'] = 0
macrophages_data_wound[is.na(macrophages_data_wound$pval),'pval'] = 1
macrophages_data_wound[is.na(macrophages_data_wound$qval),'qval'] = 1
macrophages_data_wound$feature <- rownames(macrophages_data_wound)
macrophages_data_wound <- within(macrophages_data_wound,
                                 feature <- factor(feature,
                                                   levels=order_sig_1))


stromal_data_wound <- stromal_data_wound[order_sig,]
rownames(stromal_data_wound) <- order_sig
stromal_data_wound[is.na(stromal_data_wound$coef),'coef'] = 0
stromal_data_wound[is.na(stromal_data_wound$pval),'pval'] = 1
stromal_data_wound[is.na(stromal_data_wound$qval),'qval'] = 1
stromal_data_wound$feature <- rownames(stromal_data_wound)
stromal_data_wound <- within(stromal_data_wound,
                             feature <- factor(feature,
                                               levels=order_sig_1))

fibroblast_data_wound_temp_diff_bar <- diff_bar_plot(fibroblast_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
macrophages_data_wound_temp_diff_bar <- diff_bar_plot(macrophages_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')
stromal_data_wound_temp_diff_bar <- diff_bar_plot(stromal_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')
box_association <- readRDS("matched/Tweedieverse_output_Epithelial/figures/Sample_type_gg_associations.RDS")
box_association_1 <- readRDS("matched/Tweedieverse_output_Endothelial/figures/Sample_type_gg_associations.RDS")
box_association_2 <- readRDS("matched/Tweedieverse_output_Stem/figures/Sample_type_gg_associations.RDS")


fig_3 <- ggdraw() +
  draw_plot(fibroblast_data_wound_temp_diff_bar,
            x = 0, y = .47, width = .40, height = .55) +
  draw_plot(macrophages_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .40, y = .47, width = .35, height = .55) +
  draw_plot(stromal_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .75, y = .47, width = .25, height = .55)  +
  draw_plot(box_association[[1]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association_1[[10]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association_2[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a", "Epithelial", "Endothelial" ,"Stem" , "b", "c", "d", "e")),
                  size = 7,x = c(0, .20, .40, .75, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))

ggsave(filename = 'unmatched_barplot_fibroblast.pdf', plot=fig_4, width = 183, height = 110, units = "mm", dpi = 350)

## read for fibroblast
fibroblast <- read.delim(
  "unmatched/Tweedieverse_output_Fibroblasts/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
fibroblast_data_wound <- fibroblast[fibroblast$metadata=="Sample_type" & fibroblast$value=="Wound" ,]
rownames(fibroblast_data_wound) <- fibroblast_data_wound$feature


## read for macrophges
macrophages <- read.delim(
  "unmatched/Tweedieverse_output_Macrophages/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

macrophages_data_wound <- macrophages[macrophages$metadata=="Sample_type" & macrophages$value=="Wound" ,]
rownames(macrophages_data_wound) <- macrophages_data_wound$feature


## read for stromal
stromal <- read.delim(
  "unmatched/Tweedieverse_output_Stromal/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

stromal_data_wound <- stromal[stromal$metadata=="Sample_type" & stromal$value=="Wound" ,]
rownames(stromal_data_wound) <- stromal_data_wound$feature

##Using the defined genes which are markers for different cell types###
imp_genes <- c("CD68", "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1","Krt14", "Krt10", "Krt5","Plin1", "Adipoq", "Pparg", "Fabp4","Ptprc","Pecam1", "CD34")


order_sig <- rownames(fibroblast_data_wound)[1:number_of_sig_to_keep]
order_sig <- unique(append(order_sig,imp_genes))#Put a unique before append
##add the association for each celltype
fibroblast_data_wound <- fibroblast_data_wound[order_sig,]
##reinitalizing the rownames to the list 
rownames(fibroblast_data_wound) <- order_sig
fibroblast_data_wound[is.na(fibroblast_data_wound$coef),'coef'] = 0
fibroblast_data_wound[is.na(fibroblast_data_wound$pval),'pval'] = 1
fibroblast_data_wound[is.na(fibroblast_data_wound$qval),'qval'] = 1
fibroblast_data_wound$feature <- rownames(fibroblast_data_wound)
fibroblast_data_wound<- fibroblast_data_wound[order(fibroblast_data_wound$coef),]
order_sig_1 <- rownames(fibroblast_data_wound)
fibroblast_data_wound <- within(fibroblast_data_wound,
                                feature <- factor(feature,
                                                  levels=order_sig_1))

macrophages_data_wound <- macrophages_data_wound[order_sig,]
rownames(macrophages_data_wound) <- order_sig
macrophages_data_wound[is.na(macrophages_data_wound$coef),'coef'] = 0
macrophages_data_wound[is.na(macrophages_data_wound$pval),'pval'] = 1
macrophages_data_wound[is.na(macrophages_data_wound$qval),'qval'] = 1
macrophages_data_wound$feature <- rownames(macrophages_data_wound)
macrophages_data_wound <- within(macrophages_data_wound,
                                 feature <- factor(feature,
                                                   levels=order_sig_1))


stromal_data_wound <- stromal_data_wound[order_sig,]
rownames(stromal_data_wound) <- order_sig
stromal_data_wound[is.na(stromal_data_wound$coef),'coef'] = 0
stromal_data_wound[is.na(stromal_data_wound$pval),'pval'] = 1
stromal_data_wound[is.na(stromal_data_wound$qval),'qval'] = 1
stromal_data_wound$feature <- rownames(stromal_data_wound)
stromal_data_wound <- within(stromal_data_wound,
                             feature <- factor(feature,
                                               levels=order_sig_1))

fibroblast_data_wound_temp_diff_bar <- diff_bar_plot(fibroblast_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
macrophages_data_wound_temp_diff_bar <- diff_bar_plot(macrophages_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')
stromal_data_wound_temp_diff_bar <- diff_bar_plot(stromal_data_wound, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="qval", orderby = NA, x_label = 'Effect size', y_label = '')

box_association <- readRDS("unmatched/Tweedieverse_output_Fibroblasts/figures/Sample_type_gg_associations.RDS")
box_association_1 <- readRDS("unmatched/Tweedieverse_output_Macrophages/figures/Sample_type_gg_associations.RDS")
box_association_2 <- readRDS("unmatched/Tweedieverse_output_Stromal/figures/Sample_type_gg_associations.RDS")



fig_4 <- ggdraw() +
  draw_plot(fibroblast_data_wound_temp_diff_bar,
            x = 0, y = .47, width = .40, height = .55) +
  draw_plot(macrophages_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .40, y = .47, width = .35, height = .55) +
  draw_plot(stromal_data_wound_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .75, y = .47, width = .25, height = .55)  +
  draw_plot(box_association[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[14]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association_1[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association_2[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a",  "Fibroblast", "Macrophages","Stromal" , "b", "c", "d", "e")),
                  size = 7,x = c(0, .07, .40, .77, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))

ggsave(filename = 'unmatched_barplot.pdf', plot=fig_3, width = 183, height = 110, units = "mm", dpi = 350)
