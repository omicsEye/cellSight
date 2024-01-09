#' Title
#'
#' @param int_seur
#' @param resolution
#' @param cluster_name
#'
#' @return
#' @export
#'
#' @examples
pca_clustering<-function(int_seur, output_directory, resolution = "integrated_snn_res.0.8",cluster_name= NULL){
  #setwd(output_directory)

  int_seur <- int_seur |>
    FindVariableFeatures() |>
    RunPCA()
  path <- paste0(output_directory,"/sample_plot/")
  # Check if the directory exists
  if (!dir.exists(path)) {
    # If it doesn't exist, create it
    dir.create(path, recursive = TRUE)
    cat("Directory created:", path, "\n")
  } else {
    cat("Directory already exists:", path, "\n")
  }
  file_path_1 <- paste0(path,"dimplot_type.png")
  dim_type <- DimPlot(int_seur, group.by = "type")
  ggsave(file_path_1,dim_type,bg="white")

  file_path_2 <- paste0(path,"dimplot_sample.png")
  dim_sample <- DimPlot(int_seur, group.by = "sample")
  ggsave(file_path_2,dim_sample,bg="white")

  file_path_3 <- paste0(path,"elbow_plot.png")
  elbow_plot <- ElbowPlot(int_seur, ndims = 50)
  ggsave(file_path_3, elbow_plot,bg="white")

  DefaultAssay(int_seur) <- "integrated"

  int_seur <- int_seur |>
    RunUMAP(dims = 1:30) |>
    FindNeighbors(dims = 1:30) |>
    FindClusters()
  directory_path <- paste0(output_directory,"/pca_clusters/")

  # Check if the directory exists
  if (!dir.exists(directory_path)) {
    # If it doesn't exist, create it
    dir.create(directory_path, recursive = TRUE)
    cat("Directory created:", directory_path, "\n")
  } else {
    cat("Directory already exists:", directory_path, "\n")
  }
  file<- paste0(directory_path,"/dimplots/")

  if (!dir.exists(file)) {
    # If it doesn't exist, create it
    dir.create(file, recursive = TRUE)
    cat("Directory created:", file, "\n")
  } else {
    cat("Directory already exists:", file, "\n")
  }

  dim_plot_0.2 <- DimPlot(int_seur,
          group.by = "integrated_snn_res.0.2",
          label = T,
          label.box = T) +
    theme(legend.position = "none")

  file_path_4 <- paste0(file,"integrated_snn_res(0.2).png")
  ggsave(file_path_4,dim_plot_0.2,bg="white")

  dim_plot_0.4 <- DimPlot(int_seur,
          group.by = "integrated_snn_res.0.4",
          label = T,
          label.box = T) +
    theme(legend.position = "none")

  file_path_5 <- paste0(file,"integrated_snn_res(0.4).png")
  ggsave(file_path_4,dim_plot_0.4,bg="white")

  dim_plot_0.6 <- DimPlot(int_seur,
          group.by = "integrated_snn_res.0.6",
          label = T,
          label.box = T) +
    theme(legend.position = "none")

  file_path_6 <- paste0(file,"integrated_snn_res(0.6).png")
  ggsave(file_path_6,dim_plot_0.6,bg="white")

  int_seur |>
    DimPlot(split.by = "sample")
  Idents(int_seur) <- resolution
  int_seur <- int_seur |>
    NormalizeData(assay = "RNA")
  if(!is.null(cluster_name)){
    new.cluster.ids <- cluster_name
    names(new.cluster.ids) <- levels(int_seur)
    int_seur <- RenameIdents(int_seur, new.cluster.ids)
    dim_plot_0.8 <- DimPlot(int_seur, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

  }
  else{
  dim_plot_0.8 <- DimPlot(int_seur, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  }

  file_path_7 <- paste0(file,"integrated_snn_res(0.8).png")
  ggsave(file_path_7,dim_plot_0.8,bg="white")




  for (i in int_seur$integrated_snn_res.0.8 |> unique()) {
    marker_file<- paste0(directory_path,"/markers/")
    if (!dir.exists(file)) {
      # If it doesn't exist, create it
      dir.create(file, recursive = TRUE)
      cat("Directory created:", file, "\n")
      int_seur |>
        FindConservedMarkers(ident.1 = i,
                             grouping.var = "type") |>
        write.csv(
          paste0(
            "~/analysis/all-",
            i,
            ".csv")
          )
    } else {
      cat("Directory already exists:", file, "\n")
    }

    int_seur |>
      FindConservedMarkers(ident.1 = i,
                           grouping.var = "type") |>
      write.csv(
        paste0(
          "~/analysis/all-",
          i,
          ".csv")
      )

  }
}
