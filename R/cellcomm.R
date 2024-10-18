#' Title
#'
#' @param obj_list
#' @param output_directory
#' @param imp_var
#' @param species "human" or "mouse"
#' @return
#' @import CellChat
#'
#' @export
#'
#' @examples
cellcomm_analysis <- function(obj_list, output_directory, imp_var, species = "human") {

  # If imp_var is not "all", perform analysis on a subset
  if (imp_var != "all") {
    obj_list <- pca_clusters
    imp <- subset(x = obj_list, subset = (sample == imp_var))
    data.input <- GetAssayData(imp, assay = "RNA", slot = "data")

    labels <- Idents(imp)
    meta <- data.frame(group = labels, row.names = names(labels))

    # Create CellChat object
    cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")

    # Load appropriate CellChat database based on species
    if (species == "human") {
      data("CellChatDB.human")  # Load the human CellChat database
      CellChatDB <- CellChatDB.human
      showDatabaseCategory(CellChatDB)
    } else if (species == "mouse") {
      data("CellChatDB.mouse")  # Load the mouse CellChat database
      CellChatDB <- CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
    } else {
      stop("Database not available for the specified species.")
    }

    # Set the CellChatDB for further analysis
    cellchat_normal@DB <- CellChatDB
    cellchat_normal <- subsetData(cellchat_normal)
    cellchat_normal <- identifyOverExpressedGenes(cellchat_normal, thresh.p = 1, do.fast = FALSE)
    cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

    # Clean memory
    gc()

    # Compute communication probability
    cellchat_normal <- computeCommunProb(cellchat_normal, trim = 0, nboot = 1)
    gc()

    # Filter out cell-cell communication with low cell counts
    cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)

    # Compute communication probability pathways
    cellchat_normal <- computeCommunProbPathway(cellchat_normal)
    df.netP <- reshape2::melt(cellchat_normal@net$prob, value.name = "prob")

    df_significant_netP <- df.netP[apply(df.netP != 0, 1, all),]

    # Aggregate network and set pathways
    cellchat_normal <- aggregateNet(cellchat_normal, thresh = 0.05)
    cellchat_normal@net$pathways <- df_significant_netP$Var3
    groupSize <- as.numeric(table(cellchat_normal@idents))

    # Define the directory and check if it exists
    fig_dir <- paste0(output_directory, "/cellcommunication/")

    # If the directory doesn't exist, create it
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }

    # Open the SVG device for the first plot (interaction counts)
    svglite(file = paste0(fig_dir, "pathway_count.svg"), width = 11, height = 10)

    # Set up for two plots side by side
    par(mfrow = c(1, 2), xpd = TRUE, cex.main = 1.5, font.main = 2)

    # Plot the number of interactions
    netVisual_circle(cellchat_normal@net$count,
                     vertex.weight = groupSize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = "Number of interactions",
                     vertex.label.cex = 0.5)

    # Close the device after saving the first plot
    dev.off()

    # Open the SVG device for the second plot (interaction weights/strength)
    svglite(file = paste0(fig_dir, "pathway_weight.svg"), width = 11, height = 7.2)

    # Set up for two plots side by side again
    par(mfrow = c(1, 2), xpd = TRUE, cex.main = 1.5, font.main = 2)

    # Plot the interaction weights/strength
    netVisual_circle(cellchat_normal@net$weight,
                     vertex.weight = groupSize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = "Interaction weights/strength",
                     vertex.label.cex = 0.5)

    # Close the device after saving the second plot
    dev.off()

  } else {
    # Handle case when `imp_var == "all"` (analysis for entire object)
    obj_list <- pca_clusters
    imp <- obj_list
    data.input <- GetAssayData(imp, assay = "RNA", slot = "data")

    labels <- Idents(imp)
    meta <- data.frame(group = labels, row.names = names(labels))

    # Create CellChat object
    cellchat_normal <- createCellChat(object = data.input, meta = meta, group.by = "group")

    # Load appropriate CellChat database based on species
    if (species == "human") {
      data("CellChatDB.human")  # Load the human CellChat database
      CellChatDB <- CellChatDB.human
      showDatabaseCategory(CellChatDB)
    } else if (species == "mouse") {
      data("CellChatDB.mouse")  # Load the mouse CellChat database
      CellChatDB <- CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
    } else {
      stop("Database not available for the specified species.")
    }

    # Set the CellChatDB for further analysis
    cellchat_normal@DB <- CellChatDB
    cellchat_normal <- subsetData(cellchat_normal)
    cellchat_normal <- identifyOverExpressedGenes(cellchat_normal, thresh.p = 1, do.fast = FALSE)
    cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

    # Clean memory
    gc()

    # Compute communication probability
    cellchat_normal <- computeCommunProb(cellchat_normal, trim = 0, nboot = 1)
    gc()

    # Filter out cell-cell communication with low cell counts
    cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)

    # Compute communication probability pathways
    cellchat_normal <- computeCommunProbPathway(cellchat_normal)
    df.netP <- reshape2::melt(cellchat_normal@net$prob, value.name = "prob")

    df_significant_netP <- df.netP[apply(df.netP != 0, 1, all),]

    # Aggregate network and set pathways
    cellchat_normal <- aggregateNet(cellchat_normal, thresh = 0.05)
    cellchat_normal@net$pathways <- df_significant_netP$Var3
    groupSize <- as.numeric(table(cellchat_normal@idents))

    # Define the directory and check if it exists
    fig_dir <- paste0(output_directory, "/cellcommunication/")

    # If the directory doesn't exist, create it
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }

    # Open the SVG device for the first plot (interaction counts)
    svglite(file = paste0(fig_dir, "pathway_count.svg"), width = 11, height = 10)

    # Set up for two plots side by side
    par(mfrow = c(1, 2), xpd = TRUE, cex.main = 1.5, font.main = 2)

    # Plot the number of interactions
    netVisual_circle(cellchat_normal@net$count,
                     vertex.weight = groupSize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = "Number of interactions",
                     vertex.label.cex = 0.5)

    # Close the device after saving the first plot
    dev.off()

    # Open the SVG device for the second plot (interaction weights/strength)
    svglite(file = paste0(fig_dir, "pathway_weight.svg"), width = 11, height = 7.2)

    # Set up for two plots side by side again
    par(mfrow = c(1, 2), xpd = TRUE, cex.main = 1.5, font.main = 2)

    # Plot the interaction weights/strength
    netVisual_circle(cellchat_normal@net$weight,
                     vertex.weight = groupSize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = "Interaction weights/strength",
                     vertex.label.cex = 0.5)

    # Close the device after saving the second plot
    dev.off()
  }
}
