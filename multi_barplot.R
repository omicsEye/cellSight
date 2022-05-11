library(tidyr)
library(dplyr)
library(reshape2)
library(deepath)
#setting the working directory
setwd("~/Projects/")

number_of_sig_to_keep <- 20
sig_threshold <- 0.05

## read metabolites
metabolites_Tweedieverse <- read.delim(
  "analysis/meatbolites_Tweedieverse/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_score_data_severe <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="Severe" ,]
rownames(metabolites_score_data_severe) <- metabolites_score_data_severe$feature

metabolites_score_data_non_severe <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="non-Severe" ,]
rownames(metabolites_score_data_non_severe) <- metabolites_score_data_non_severe$feature

metabolites_score_data_non_covid <- metabolites_Tweedieverse[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="non-COVID-19" ,]
rownames(metabolites_score_data_non_covid) <- metabolites_score_data_non_covid$feature



# use score_data_severe is reference
order_sig <- rownames(metabolites_score_data_severe)[1:number_of_sig_to_keep]
metabolites_score_data_severe <- metabolites_score_data_severe[order_sig,]
metabolites_score_data_severe<- metabolites_score_data_severe[order(metabolites_score_data_severe$coef),]
order_sig <- rownames(metabolites_score_data_severe)
metabolites_score_data_severe <- within(metabolites_score_data_severe,
                                        feature <- factor(feature,
                                                          levels=order_sig))
metabolites_score_data_non_severe <- metabolites_score_data_non_severe[rownames(metabolites_score_data_severe),]
metabolites_score_data_non_severe <- within(metabolites_score_data_non_severe,
                                            feature <- factor(feature,
                                                              levels=order_sig))


metabolites_score_data_non_covid <- metabolites_score_data_non_covid[rownames(metabolites_score_data_severe),]
metabolites_score_data_non_covid <- within(metabolites_score_data_non_covid,
                                           feature <- factor(feature,
                                                             levels=order_sig))

metabolites_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_severe, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
metabolites_non_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_severe, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
metabolites_non_covid_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_covid, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')



## read association
box_association <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-Omics/analysis/meatbolites_Tweedieverse/figures/Group_gg_associations.RDS")
## do plots

fig2_metabolites <- ggdraw() +
  draw_plot(metabolites_severe_temp_diff_bar,
            x = 0, y = .47, width = .55, height = .53) +
  draw_plot(metabolites_non_severe_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .55, y = .47, width = .225, height = .53) +
  draw_plot(metabolites_non_covid_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.line.y = element_blank()),
            x = .775, y = .47, width = .225, height = .53) +
  draw_plot(box_association[[11]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[52]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[139]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[2]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a",  "Severe", "non-Severe", "non-COVID", "b", "c", "d", "e")),
                  size = 7,x = c(0, .28, .53, .76, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))
fig3_metabolites

ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'figures/fig3/fig#_barplot.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)


