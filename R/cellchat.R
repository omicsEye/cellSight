devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sqjin/CellChat")
library(CellChat)
data.input <- GetAssayData(pos_cells_nonwound_macro_pdgfrb, assay = "RNA", slot = "data") 
Idents(pos_cells_nonwound_macro_pdgfrb) <- pos_cells_nonwound_macro_pdgfrb@active.ident
new.cluster.ids <- c("1","2","3","4","5","6","7","8")
names(new.cluster.ids) <- levels(pos_cells_nonwound_macro_pdgfrb)
pos_cells_nonwound_macro_pdgfrb <- RenameIdents(pos_cells_nonwound_macro_pdgfrb, new.cluster.ids)
labels <- Idents(pos_cells_nonwound_macro_pdgfrb)
meta <- data.frame(group = labels, row.names = names(labels)) 
# create a CellChat object
cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")

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

cellchat_normal <- computeCommunProb(cellchat_normal)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)


cellchat_normal <- computeCommunProbPathway(cellchat_normal)
cellchat_normal <- aggregateNet(cellchat_normal)

groupSize <- as.numeric(table(cellchat_normal@idents))
par(mfrow = c(1,2), xpd=TRUE)
count <- netVisual_circle(cellchat_normal@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

ggsave(filename = 'count_cellchat.pdf', plot=count, width = 181, height = 110, units = "mm", dpi = 350)

weight <- netVisual_circle(cellchat_normal@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



mat <- cellchat_normal@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
