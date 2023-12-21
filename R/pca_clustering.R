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
pca_clustering<-function(int_seur, resolution = "integrated_snn_res.0.8",cluster_name= NULL){
    DefaultAssay(int_seur) <- "integrated"

  int_seur <- int_seur |>
    FindVariableFeatures() |>
    RunPCA()|>
    RunUMAP(dims = 1:30) |>
    FindNeighbors(dims = 1:30) |>
    FindClusters()
  Idents(int_seur) <- resolution
  int_seur <- int_seur |>
    NormalizeData(assay = "RNA")
  if(!is.null(cluster_name)){
    new.cluster.ids <- cluster_name
    names(new.cluster.ids) <- levels(int_seur)
    int_seur <- RenameIdents(int_seur, new.cluster.ids)
  }
  for (i in int_seur$integrated_snn_res.0.8 |> unique()) {
    if (!dir.exists("~/analysis/")){
      dir.create("~/analysis/")
      int_seur |>
        FindConservedMarkers(ident.1 = i) |>
        write.csv(
          paste0(
            "~/analysis/conserved-",
            i,
            ".csv"
          )
        )
    }else{
      int_seur |>
        FindConservedMarkers(ident.1 = i) |>
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
