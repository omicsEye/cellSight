
#' Quality control filtering with plots
#'
#' @param seur_object The seurat object from the previous directory step
#'
#' @return QC plots
#' @export
#'
#' @examples

qc_plots<- function(seur_object,output_directory){
  setwd(output_directory)
  obj_list <- seur_object
  for(index in 1:length(obj_list)){
    print(index)
    obj_list[[index]][["percent.mt"]] <- obj_list[[index]] %>%
      PercentageFeatureSet(pattern = "^mt-")

    plot_violin <- VlnPlot(obj_list[[index]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                           ncol = 3)
    file_path <- paste0(output_directory,"qcplot_violin_",index,".png")
    ggsave(file_path,plot_violin,bg="white")
    plot_scatter <- FeatureScatter(obj_list[[index]],feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    file_path_1 <- paste0(output_directory,"qcplot_scatter_",index,".png")
    ggsave(file_path_1,plot_scatter,bg="white")
    grid <- plot_grid(
      plot_violin, plot_scatter,
      labels = "AUTO", ncol = 1
    )
    file_path_2 <- paste0(output_directory,"qcplot_grid_joined_",index,".png")
    save_plot(file_path_2,grid,bg="white")
  }


}
