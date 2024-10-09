#' Title
#'
#' @param obj_list
#'
#' @return
#' @export
#'
#' @examples
cellcomm_analysis<-function(obj_list,output_directory,imp_var,species){
  if(imp_var != "all"){
    obj_list <- pca_clusters
    imp <- subset(x = obj_list, subset = (sample  == imp_var))
    data.input <- GetAssayData(imp, assay = "RNA", slot = "data")

    labels <- Idents(imp)
    meta <- data.frame(group = labels, row.names = names(labels))

    cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")
    if(species == "human")
    {
      CellChatDB <- CellChatDB.human
      showDatabaseCategory(CellChatDB)
    }
    if(species == "mouse")
    {
      CellChatDB <- CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
    }
    else{
      cat("Database not available")
    }
    CellChatDB.use <- CellChatDB
    cellchat_normal@DB <- CellChatDB.use
    cellchat_normal <- subsetData(cellchat_normal)
    cellchat_normal <- identifyOverExpressedGenes(cellchat_normal,thresh.p = 1,do.fast = F)
    cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

    gc()

    cellchat_normal <- computeCommunProb(cellchat_normal,trim = 0, nboot = 1)

    gc()

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)


    cellchat_normal <- computeCommunProbPathway(cellchat_normal)
    df.netP <- reshape2::melt(cellchat_normal@net$prob, value.name = "prob")

    df_significant_netP <- df.netP[apply(df.netP!=0, 1, all),]


    cellchat_normal <- aggregateNet(cellchat_normal,thresh = 0.05)
    cellchat_normal@net$pathways <- df_significant_netP$Var3
    groupSize <- as.numeric(table(cellchat_normal@idents))
  }
  else{
    obj_list <- pca_clusters
    imp <- obj_list
    data.input <- GetAssayData(imp, assay = "RNA", slot = "data")

    labels <- Idents(imp)
    meta <- data.frame(group = labels, row.names = names(labels))

    cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")
    if(species == "human")
    {
      CellChatDB <- CellChatDB.human
      showDatabaseCategory(CellChatDB)
    }
    if(species == "mouse")
    {
      CellChatDB <- CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
    }
    else{
      cat("Database not available")
    }
    CellChatDB.use <- CellChatDB
    cellchat_normal@DB <- CellChatDB.use
    cellchat_normal <- subsetData(cellchat_normal)
    cellchat_normal <- identifyOverExpressedGenes(cellchat_normal,thresh.p = 1,do.fast = F)
    cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

    gc()

    cellchat_normal <- computeCommunProb(cellchat_normal,trim = 0, nboot = 1)

    gc()

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)


    cellchat_normal <- computeCommunProbPathway(cellchat_normal)
    df.netP <- reshape2::melt(cellchat_normal@net$prob, value.name = "prob")

    df_significant_netP <- df.netP[apply(df.netP!=0, 1, all),]


    cellchat_normal <- aggregateNet(cellchat_normal,thresh = 0.05)
    cellchat_normal@net$pathways <- df_significant_netP$Var3
    groupSize <- as.numeric(table(cellchat_normal@idents))
  }
}
