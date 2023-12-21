#' Title
#'
#' @param obj_list
#'
#' @return
#' @export
#'
#' @examples
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
      output = paste0('~/tweedieverse_run/cluster_', i),
      prev_threshold = 0.0,
      entropy_threshold = 0.0,
      base_model = 'CPLM',
      plot_heatmap = T,
      plot_scatter = T,
      reference = c("type, Normal")
    )
  }
}
