
#' Title
#'
#' @param obj_list
#'
#' @return
#' @export
#'
#' @examples

sctransform_integration<-function(obj_list,output_directory){
  #setwd(output_directory)
  if(length(obj_list) == 1 ){
    obj_list <- obj_list %>%
      lapply(FUN = function(x) {
        SCTransform(x, vst.flavor = "v2")
      })
    directory_path <- paste0(output_directory,"/sctransform.rds")
    if (!dir.exists(directory_path)) {
      # If it doesn't exist, create it
      dir.create(directory_path, recursive = TRUE)
      cat("Directory created:", directory_path, "\n")
    } else {
      cat("Directory already exists:", directory_path, "\n")
    }
    #file <- paste0(output_directory,"sctransform.rds")
    saveRDS(obj_list,directory_path)
    return(obj_list)
  }
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

    directory_path <- paste0(output_directory,"/sctransform/")
    if (!dir.exists(directory_path)) {
      # If it doesn't exist, create it
      dir.create(directory_path, recursive = TRUE)
      cat("Directory created:", directory_path, "\n")
    } else {
      cat("Directory already exists:", directory_path, "\n")
    }
    file <- paste0(output_directory,"sctransform.rds")
    saveRDS(obj_list,file)

  }
  file <- paste0(output_directory,"sctransform.rds")
  saveRDS(obj_list,file)
  return(obj_list)
}
