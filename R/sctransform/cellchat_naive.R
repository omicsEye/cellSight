library(CellChat)
library(gridExtra)
library(ggpubr)
library(svglite)
theme_set(theme_pubr())

data.input <- GetAssayData(naive, assay = "RNA", slot = "data") 

labels <- Idents(naive)
meta <- data.frame(group = labels, row.names = names(labels)) 

cellchat_naive <- createCellChat(object = data.input, meta = meta, group.by = "group")

rm(seur_obj,obj_list,meta,wounded,data.input,p,p1,p2,p3,umap_theme,cluster.markers)

##import cellchat database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat_naive@DB <- CellChatDB.use


cellchat_naive <- subsetData(cellchat_naive) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4)


cellchat_naive <- identifyOverExpressedGenes(cellchat_naive,thresh.p = 1)
cellchat_naive <- identifyOverExpressedInteractions(cellchat_naive)
#df.net <- subsetCommunication(cellchat_naive,slot.name = "net")

#feature.name = as.data.frame(rownames(wounded))
#net <- netMappingDEG(cellchat_naive, features.name = df_try)
#net.up <- subsetCommunication(cellchat_naive, net = cellchat_naive@net, ligand.logFC = 0.2, receptor.logFC = NULL)

#net.up <- subsetCommunication(cellchat, net = netP, datasets = "NPD",ligand.logFC = 0.2, receptor.logFC = NULL)

#net.down <- subsetCommunication(cellchat, net = net, datasets = "LPD",ligand.logFC = -0.1, receptor.logFC = -0.1)

gc()

cellchat_naive <- computeCommunProb(cellchat_naive,trim = 0, nboot = 1)

gc()

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_naive <- filterCommunication(cellchat_naive, min.cells = 10)


cellchat_naive <- computeCommunProbPathway(cellchat_naive)


df_significant_netP <- df.netP[apply(df.netP!=0, 1, all),]


cellchat_naive <- aggregateNet(cellchat_naive,thresh = 0.05)
cellchat_naive@net$pathways <- df_significant_netP$Var3
groupSize <- as.numeric(table(cellchat_naive@idents))
par(mfrow = c(1,2), xpd=TRUE)
count <- netVisual_circle(cellchat_naive@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


weight <- netVisual_circle(cellchat_naive@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

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


mat <- cellchat_naive@net$weight
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
pathways.show.all <- cellchat_naive@net$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_naive@idents)
vertex.receiver = c(1,2,3,5,12,17)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  gg1 <- netVisual(cellchat_naive, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg2 <- netVisual_aggregate(cellchat_naive, signaling = pathways.show.all[i], layout = "chord")
  gg <- netAnalysis_contribution(cellchat_naive, signaling = pathways.show.all[i])
  #ggsave(filename=paste0(pathways.show.all[i], "_wound2_L-R_contribution.svg"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  #ggsave(filename=paste0(pathways.show.all[i], "_net_visual.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
  #ggsave(filename=paste0(pathways.show.all[i], "_wound1_net_contribution.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
  
}

i = 1
gg <- netAnalysis_contribution(cellchat_naive, signaling = pathways.show.all[i])
netVisual(cellchat_naive, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")



pathways.show.all <- cellchat_naive@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_naive@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat_naive, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat_naive, signaling = pathways.show.all[i])
  #ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat_naive, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

#> Comparing communications on a single object


# Compute the network centrality scores
cellchat_naive <- netAnalysis_computeCentrality(cellchat_naive, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_naive, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_naive)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_naive)
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

