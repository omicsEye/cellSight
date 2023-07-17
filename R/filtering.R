#' Title
#'
#' @param obj_list
#' @param index
#' @param min_rna
#' @param nfeature_rna
#' @param ncount_rna
#' @param mit.percent
#' @param sample_type
#' @param sample_name
#' @param stepwise
#'
#' @return
#' @export
#'
#' @examples

filtering<-function(obj_list,index=1,min_rna = 200,nfeature_rna = 2500,ncount_rna = 8000,mit.percent = 1,sample_type=NULL,sample_name=NULL,stepwise = FALSE){
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
  #obj_list[[index]]$sample <- sample_name
  #obj_list[[index]]$type <- sample_type
  return(obj_list)
}
