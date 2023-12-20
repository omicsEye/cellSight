

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
#'
#' @return The output form [print()]
#' @export
#'
#' @examples
data_directory<- function(dir){
  #box_dir <- "~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data/"
  #box_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/"
  box_dir <- dir
  dirs <- list.dirs(path = box_dir,
                    recursive = F,
                    full.names = F)



  obj_list <- dirs %>%
    dplyr::set_names() %>%
    { map(.f = function(x) {
      paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
        Read10X() %>%
        CreateSeuratObject()
    })}

  ## Just changed the length to dimension to get both the size##
  #dim(obj_list$Nonwound1),dim(obj_list$Nonwound2),dim(obj_list$Wound1),dim(obj_list$Wound2)
  for(i in 1:length(obj_list)){print(dim(obj_list[[i]]))}
  return(obj_list)
}
