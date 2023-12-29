
#' Title
#'
#' @param obj_list
#'
#' @return
#' @export
#'
#' @examples

sctransform_integration<-function(obj_list,output_directory){
  setwd(output_directory)
  if(length(obj_list) == 1 ){
    obj_list <- obj_list %>%
      lapply(FUN = function(x) {
        SCTransform(x, vst.flavor = "v2")
      })
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
    return(obj_list)
  }
  saveRDS(
    "~/sc_transform/integrated.rds"
  )
}
