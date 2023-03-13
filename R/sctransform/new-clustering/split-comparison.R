library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
my_cols = brewer.pal(4,"Dark2")
int <-
  read_rds("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/objects/sc-integrated.rds")

nonint <-
  read_rds(
    "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/objects/sc-non-integrated.rds"
  )



nonint$type <-
  nonint$sample %>%
  gsub('[[:digit:]]+', '', .)

int$type <-
  int$sample %>%
  gsub('[[:digit:]]+', '', .)


nonint_sample <- nonint %>%
  DimPlot(group.by = "sample",cols=alpha(my_cols,0.4),pt.size=0.3) +
  labs(title = "Before integration") +
  theme(axis.text = element_blank()) +
   NoAxes()+ggtitle(" ")+NoLegend()+
  scale_fill_brewer(palette = "Dark2")

nonint_injury <- nonint %>%
  DimPlot(group.by = "type",cols=alpha(my_cols,0.4),pt.size=0.3) +
  labs(title = "Before integration") +
  theme(axis.text = element_blank()) +
  NoAxes()+ggtitle(" ")+NoLegend()+
  scale_fill_brewer(palette = "Dark2")

int_sample <- int %>%
  DimPlot(group.by = "sample",cols=alpha(my_cols,0.4),pt.size=0.3) +
  labs(title = "After integration") +
  theme(axis.text = element_blank()) +
  NoAxes()+ggtitle(" ")+guides(color = guide_legend(override.aes = list(size=4), ncol=1) )+
  scale_fill_brewer(palette = "Dark2")

int_injury <- int %>%
  DimPlot(group.by = "type") +
  labs(title = "After integration",cols=alpha(my_cols,0.4),pt.size=0.3) +
  theme(axis.text = element_blank()) +
  NoAxes()+ggtitle(" ")+ guides(color = guide_legend(override.aes = list(size=4), ncol=1) )+
  scale_fill_brewer(palette = "Dark2")

sc_int <- seur_obj %>%
  DimPlot(
    group.by = "Celltype",
    label = T,
    label.size = 2.5
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_blank()) +
  NoAxes()
  

sc_nonint <- seur_obj %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.size = 2.5
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text = element_blank()) +
  NoAxes()

x <- (nonint_sample +int_sample ) 

y <-  (nonint_injury+int_injury  )

fig1 <- (y|x)/ sc_int 
ggsave("C:/Users/ranoj/Desktop/Single_cell_output/plot_integrated.png", width = 6,   height = 2.5)
