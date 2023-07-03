# load libraries ----------------------------------------------------------

install.packages("rlang")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")
install.packages("dplyr")  

library(dplyr)
library(glmGamPoi)
library(DelayedMatrixStats)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(cowplot)

#Please have the data set in the format data_name/outs/filtered_feature_bc_matrix 
#and the barcodes,features and matrix inside them
##for example sample_data/outs/filtered_feature_bc_matrix/barcodes.tsv
##for example sample_data/outs/filtered_feature_bc_matrix/features.tsv
##for example sample_data/outs/filtered_feature_bc_matrix/matrix.mtx

#Testing
#local_dir <- "~/Desktop/Single_Cell_output/"
data_directory<- function(dir){
  #box_dir <- "~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data/"
  #box_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/"
  box_dir <- dir
  dirs <- list.dirs(path = box_dir, 
                    recursive = F, 
                    full.names = F)
  
 
  
  obj_list <- dirs %>%
    set_names() %>%
    map(.f = function(x) {
      paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
        Read10X() %>%
        CreateSeuratObject()
    })
  
  ## Just changed the length to dimension to get both the size##
  #dim(obj_list$Nonwound1),dim(obj_list$Nonwound2),dim(obj_list$Wound1),dim(obj_list$Wound2)
  for(i in 1:length(obj_list)){print(dim(obj_list[[i]]))}
  return(obj_list)
}

# nonwound 1 qc -----------------------------------------------------------
###Enter the identity index you want the qc plot to be displayed(Generates both scatter & Violinplot)

qc_plot <- function(seur_object){
  obj_list <- seur_object
  for(index in 1:length(obj_list)){
    print(index)
    obj_list[[index]][["percent.mt"]] <- obj_list[[index]] %>%
      PercentageFeatureSet(pattern = "^mt-")
    
    plot_violin <- VlnPlot(obj_list[[index]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           ncol = 3)
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_violin_",index,".png"))
    plot_scatter <- FeatureScatter(obj_list[[index]],feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_scatter_",index,".png"),bg="white")
    grid <- plot_grid(
      plot_violin, plot_scatter,
      labels = "AUTO", ncol = 1
    )
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_grid_joined_",index,".png"),bg="white")
  }
}

# Stepwise denotes if all qc needs to be done seperately stepwise = True, if not all
# qc will be done at the same time. Default = True.
#sample_type represent 


filtering<-function(object_list,index=1,min_rna = 200,nfeature_rna = 2500,ncount_rna = 8000,mit.percent = 1,sample_type,sample_name,stepwise = FALSE){
  for(i in 1:length(obj_list)){
    obj_list[[i]][["percent.mt"]] <- obj_list[[i]] %>%
      PercentageFeatureSet(pattern = "^mt-") 
  }
  
  if(stepwise == FALSE)  {
    for(index in 1: length(obj_list)){
      obj_list[[index]] <- 
        obj_list[[index]] %>%
        subset(nFeature_RNA > min_rna & 
                 nFeature_RNA < nfeature_rna & 
                 nCount_RNA < ncount_rna &
                 percent.mt < mit.percent)
    }
    
  }
  if(stepwise == TRUE){
    obj_list[[index]] <- 
      obj_list[[index]] %>%
      subset(nFeature_RNA > min_rna & 
               nFeature_RNA < nfeature_rna & 
               nCount_RNA < ncount_rna &
               percent.mt < mit.percent)
    
  }
  obj_list[[index]]$sample <- sample_name 
  obj_list[[index]]$type <- sample_type           
  return(obj_list)
}

sctransform_integration<-function(obj_list){
  if(length(obj_list) == 1 ){
    obj_list <- obj_list %>%
      lapply(FUN = function(x) {
        SCTransform(x, vst.flavor = "v2")
    })
  }  
  return(obj_list)
  if(length(obj_list)>1){
    obj_list <- obj_list %>%
      lapply(FUN = function(x) {
        SCTransform(x, vst.flavor = "v2")
      })
    features <- obj_list %>%
      SelectIntegrationFeatures(nfeatures = 3000)
    obj_list <- obj_list %>%
      PrepSCTIntegration(anchor.features = features)
    anchors <- obj_list %>% # ~20 min with m1
      FindIntegrationAnchors(normalization.method = "SCT",
                             anchor.features = features)
    
    obj_list <- anchors %>%
      IntegrateData(normalization.method = "SCT")
  }
  saveRDS(
    "/Users/ranoj/Desktop/Single_cell_output/integrated.rds"
  )
  return(obj_list)
}

pca_clustering<-function(obj_list, resolution = "integrated_snn_res.0.8",cluster_name= NULL){
  if (length(obj_list) == 1){
    DefaultAssay(seur_obj) <- "SCT"
  }
  else{
    DefaultAssay(seur_obj) <- "integrated"
  }
  seur_obj <- seur_obj |>
    RunUMAP(dims = 1:30) |>
    FindNeighbors(dims = 1:30) |>
    FindClusters()
  Idents(seur_obj) <- resolution
  seur_obj <- seur_obj |>
    NormalizeData(assay = "RNA")
  if(!is.null(cluster_name)){
    new.cluster.ids <- cluster_name
    names(new.cluster.ids) <- levels(seur_obj)
    seur_obj <- RenameIdents(seur_obj, new.cluster.ids)
  }
  for (i in seur_obj$integrated_snn_res.0.8 |> unique()) {
    if (!dir.exists("~/analysis/")){
      dir.create("~/analysis/")
      seur_obj |>
        FindConservedMarkers(ident.1 = i,
                             grouping.var = "type") |>
        write.csv(
          paste0(
            "~/analysis/conserved-",
            i,
            ".csv"
          )
        )
    }else{
      seur_obj |>
        FindConservedMarkers(ident.1 = i,
                             grouping.var = "type") |>
        write.csv(
          paste0(
            "~/analysis/conserved-",
            i,
            ".csv"
          )
        )
    }
    
  }
}
tweedieverse_analysis<-function(obj_list){
  for (i in obj_list$Celltype %>% unique()) {
    print(i)
    obj_sub <- obj_list %>%
      subset(Celltype == i)
    
    counts <- obj_sub %>%
      GetAssayData(assay = "RNA", slot = "counts")
    
    input_features <- counts |>
      as.matrix() |>
      t() |>
      as.data.frame()
    
    test <- Tweedieverse(
      input_features,
      obj_list@meta.data[6] ,
      output = paste0('C:/Users/ranoj/Desktop/Single_cell_output/tweedieverse_rerun/cluster_', i),
      prev_threshold = 0.0,
      entropy_threshold = 0.0,
      base_model = 'CPLM',
      plot_heatmap = T,
      plot_scatter = T,
      reference = c("type, Normal")
    )
  }

cellchat_eval <- function(obj_list,organism){
  data.input <- GetAssayData(wounded, assay = "RNA", slot = "data") 
  
  labels <- Idents(wounded)
  meta <- data.frame(group = labels, row.names = names(labels)) 
  
  cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  rm(seur_obj,obj_list,meta,wounded,data.input,p,p1,p2,p3,umap_theme,cluster.markers)
  
  ##import cellchat database
  if(organism == "Mouse"){
  CellChatDB <- CellChatDB.mouse
  showDatabaseCategory(CellChatDB)
  
  CellChatDB.use <- CellChatDB 
  }
  if(organism == "Human"){
    CellChatDB <- CellChatDB.human
    showDatabaseCategory(CellChatDB)
    
    CellChatDB.use <- CellChatDB 
  }  
  # set the used database in the object
  cellchat_normal@DB <- CellChatDB.use
  
  
  cellchat_normal <- subsetData(cellchat_normal) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4)
  
  
  cellchat_normal <- identifyOverExpressedGenes(cellchat_normal,thresh.p = 1)
  cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)
  #df.net <- subsetCommunication(cellchat_normal,slot.name = "net")
  
  #feature.name = as.data.frame(rownames(wounded))
  #net <- netMappingDEG(cellchat_normal, features.name = df_try)
  #net.up <- subsetCommunication(cellchat_normal, net = cellchat_normal@net, ligand.logFC = 0.2, receptor.logFC = NULL)
  
  #net.up <- subsetCommunication(cellchat, net = netP, datasets = "NPD",ligand.logFC = 0.2, receptor.logFC = NULL)
  
  #net.down <- subsetCommunication(cellchat, net = net, datasets = "LPD",ligand.logFC = -0.1, receptor.logFC = -0.1)
  
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
  pathways.show.all <- cellchat_normal@net$pathways
  # check the order of cell identity to set suitable vertex.receiver
  levels(cellchat_normal@idents)
  vertex.receiver = c(1,2,3,5,12,17)
  for (i in 1:length(pathways.show.all)) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    gg1 <- netVisual(cellchat_normal, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    gg2 <- netVisual_aggregate(cellchat_normal, signaling = pathways.show.all[i], layout = "chord")
    gg <- netAnalysis_contribution(cellchat_normal, signaling = pathways.show.all[i])
    #ggsave(filename=paste0(pathways.show.all[i], "_wound2_L-R_contribution.svg"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
    #ggsave(filename=paste0(pathways.show.all[i], "_net_visual.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
    #ggsave(filename=paste0(pathways.show.all[i], "_wound1_net_contribution.svg"), plot=gg1, width = 3, height = 2, units = 'in', dpi = 300)
    
  }
  
  i = 1
  gg <- netAnalysis_contribution(cellchat_normal, signaling = pathways.show.all[i])
  netVisual(cellchat_normal, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  
  
  # Access all the signaling pathways showing significant communications
  pathways.show.all <- cellchat_normal@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  levels(cellchat_normal@idents)
  vertex.receiver = seq(1,4)
  for (i in 1:length(pathways.show.all)) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    netVisual(cellchat_normal, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    gg <- netAnalysis_contribution(cellchat_normal, signaling = pathways.show.all[i])
    #ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  }
  
  
  # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  netVisual_bubble(cellchat_normal, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
  
  #> Comparing communications on a single object
  
  
  # Compute the network centrality scores
  cellchat_normal <- netAnalysis_computeCentrality(cellchat_normal, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  netAnalysis_signalingRole_network(cellchat_normal, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat_normal)
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat_normal)
  #> Signaling role analysis on the cell-cell communication network from user's input
  gg1 + gg2
}  
}

cellsight<- function(dir){
  directory = readline(prompt = "Enter the path to your directory : ")
  seur_obj <- data_directory(directory)
  qc_plot(seur_obj)
  var1 = readline(prompt = "Enter the minimun RNA : ")
  var2 = readline(prompt = "Enter the number of feature : ")
  var3 = readline(prompt = "Enter the number of count : ")
  var4 =  readline(prompt  = "Enter the mitochondrial percentage : ")
  seur_obj <- filtering(seur_obj,min_rna = var1,nfeature_rna = var2,ncount_rna = var3,mit.percent = var4)
  int_seur <- sctransform_integration(seur_obj)
  int_seur <- pca_clustering(int_seur)
  tweedieverse_analysis(int_seur)
  cellchat_eval(int_seur,organism)
  
}



