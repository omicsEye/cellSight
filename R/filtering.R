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
#' @import Tweedieverse
#' @import Seurat
#' @import CellChat
#' @examples
#' filtering(seurat_object,output_directory)
#' filtering(seur_obj,"C:/Users/Desktop/output")
#' Filters each dataset to a pre-defined values
#' Returns a filtered dataset with the desired values

filtering<-function(obj_list,output_directory,index=1,min_rna = 200,nfeature_rna = 2500,ncount_rna = 8000,mit.percent = 1,sample_type=NULL,sample_name=NULL,stepwise = FALSE){
  #setwd(output_directory)
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
  directory_path <- paste0(output_directory,"/filtering/")

  # Check if the directory exists
  if (!dir.exists(directory_path)) {
    # If it doesn't exist, create it
    dir.create(directory_path, recursive = TRUE)
    cat("Directory created:", directory_path, "\n")
  } else {
    cat("Directory already exists:", directory_path, "\n")
  }
  file <- paste0(directory_path,"filtering.rds")
  saveRDS(obj_list,file)
  return(obj_list)
}
