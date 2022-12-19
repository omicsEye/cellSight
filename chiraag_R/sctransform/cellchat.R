library(CellChat)
library(gridExtra)
library(ggpubr)
library(svglite)
theme_set(theme_pubr())

data.input <- GetAssayData(wounded, assay = "RNA", slot = "data") 

labels <- Idents(wounded)
meta <- data.frame(group = labels, row.names = names(labels)) 

cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")

rm(seur_obj,obj_list,meta,wounded,data.input,p,p1,p2,p3,umap_theme,cluster.markers)

##import cellchat database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat_normal@DB <- CellChatDB.use


cellchat_normal <- subsetData(cellchat_normal) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4)


cellchat_normal <- identifyOverExpressedGenes(cellchat_normal)
cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

gc()

cellchat_normal <- computeCommunProb(cellchat_normal)

gc()

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)


cellchat_normal <- computeCommunProbPathway(cellchat_normal,thresh = 1)
cellchat_normal <- aggregateNet(cellchat_normal,thresh = 1)

groupSize <- as.numeric(table(cellchat_normal@idents))
par(mfrow = c(1,2), xpd=TRUE)
count <- netVisual_circle(cellchat_normal@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


weight <- netVisual_circle(cellchat_normal@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()
fig_dir <- "C:/Users/ranoj/Desktop/Single_cell_output/cellchat/"
svglite(paste0(fig_dir,file = "pathway_count.svg"), width=20, height=10)
#pdf(paste0(fig_dir, "count.pdf"), width=2.45, height=1.3)
#(p1|p2)/p3 +umap_theme
#CombinePlots(, ncol=3)
count
dev.off()

fig_dir <- "C:/Users/ranoj/Desktop/Single_cell_output/cellchat/"
svglite(paste0(fig_dir,file = "pathway_weight.svg"), width=20, height=10)
#pdf(paste0(fig_dir, "count.pdf"), width=2.45, height=1.3)
#(p1|p2)/p3 +umap_theme
#CombinePlots(, ncol=3)
weight
dev.off()

mat <- cellchat_normal@net$weight
#svg(paste0("comp_pmfs", ".svg"), width=5.3, height=7)
par(mar=c(0.5, 0.5, 0.2, 0.2), mfrow=c(4,3),
    oma = c(4, 4, 0.2, 0.2))

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

###Just change the output directory based on the sample##
setwd("C:/Users/ranoj/Desktop/Single_cell_output/cellchat/wound2")
pathways.show.all <- cellchat_normal@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_normal@idents)
vertex.receiver = c(1,2,3,5,12,17)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  gg1 <- netVisual(cellchat_normal, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  #gg2 <- netVisual_aggregate(cellchat_normal, signaling = pathways.show.all[i], layout = "chord")
  gg <- netAnalysis_contribution(cellchat_normal, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_wound2_L-R_contribution.svg"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  #ggsave(filename=paste0(pathways.show.all[i], "_net_visual.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
  #ggsave(filename=paste0(pathways.show.all[i], "_wound1_net_contribution.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
  
}
i = 1
gg <- netAnalysis_contribution(cellchat_normal, signaling = pathways.show.all[i])
netVisual(cellchat_normal, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")

ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
