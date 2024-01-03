

#' @title Directory function
#'
#' @description The function looks for the desired files under the supplied directory
#'
#' @param dir The directory path
#'
#' @import dplyr
#' @import glmGamPoi
#' @import DelayedMatrixStats
#' @import Seurat
#' @import ggplot2
#' @import tidyverse
#' @import gridExtra
#' @import cowplot
#' @import purrr
#' @import metap
#' @return The output form [print()]
#' @export
#'
#' @examples
data_directory<- function(dir,output_dir){
  #box_dir <- "~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data/"
  #box_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/"
  #setwd(output_dir)
  box_dir <- dir
  dirs <- list.dirs(path = box_dir,
                    recursive = F,
                    full.names = F)



  obj_list <- dirs %>%
    purrr::set_names() %>%
    { map(.x =.,.f = function(x) {
      paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
        Read10X() %>%
        CreateSeuratObject()
    })}

  ## Just changed the length to dimension to get both the size##
  #dim(obj_list$Nonwound1),dim(obj_list$Nonwound2),dim(obj_list$Wound1),dim(obj_list$Wound2)
  #for(i in 1:length(obj_list)){print(dim(obj_list[[i]]))}
  for (i in seq_along(obj_list)) {
    seuratObjectName <- names(obj_list)[i]
    print(seuratObjectName)
    print(dim(obj_list[[i]]))
    onlyCharacters <- gsub("[^a-zA-Z]", "", as.character(seuratObjectName))
    #print(onlyCharacters)
    obj_list[[seuratObjectName]]$sample <- seuratObjectName
    obj_list[[seuratObjectName]]$type <- onlyCharacters

  }

  print("Entering the directory creation section")
  # Specify the directory path
  directory_path <- paste0(output_dir,"data")

  # Check if the directory exists
  if (!dir.exists(directory_path)) {
    # If it doesn't exist, create it
    dir.create(directory_path, recursive = TRUE)
    cat("Directory created:", directory_path, "\n")
  } else {
    cat("Directory already exists:", directory_path, "\n")
  }
  file <- paste0(directory_path,"data_curation.rds")
  saveRDS(obj_list,file)
  return(obj_list)
}
