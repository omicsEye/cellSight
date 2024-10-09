#' Title
#'
#' @param obj_list
#'
#' @return
#' @export
#' @import Tweedieverse
#' @examples
de_analysis<-function(obj_list,output_directory,imp_var){
  imp_var = imp_var
  cat("Important Variable selected:", imp_var, "\n")
  path <- paste0(output_directory,"/tweedieverse/")
  # Check if the directory exists
  if (!dir.exists(path)) {
    # If it doesn't exist, create it
    dir.create(path, recursive = TRUE)
    cat("Directory created:", path, "\n")
  } else {
    cat("Directory already exists:", path, "\n")
  }
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

    test <- Tweedieverse::Tweedieverse(
      input_features,
      obj_sub@meta.data[imp_var] ,
      output = paste0(path,'/cluster_', i),
      prev_threshold = 0.0,
      entropy_threshold = 0.0,
      base_model = 'CPLM',
      plot_heatmap = T,
      plot_scatter = T,
    )
  }
}
