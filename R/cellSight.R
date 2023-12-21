#' cellSight Function
#'
#' This function performs a series of analysis steps on single-cell data using the
#' following functions in order: data_directory, filtering, pca_clusters, and qc_plots.
#'
#' @param input_data A data frame or matrix containing the single-cell data.
#' @param other_parameters Additional parameters to be passed to the sub-functions.
#'
#' @return A list containing the results of each analysis step.
#'
#' @export
#'
#' @examples
#' data(iris)
#' result <- cellSight(iris)
#'
#'
#'
#' @seealso
#' \code{\link{data_directory}}, \code{\link{filtering}}, \code{\link{pca_clusters}}, \code{\link{qc_plots}}
#'
#' @family analysis
#' @keywords single-cell data
#' @aliases cellSight
#' @concept Cell analysis
#' @concept Single-cell RNA-seq
cellSight <- function(data_directory,output_directory) {
  # Step 1: Call data_directory function
  data_dir_result <- data_directory(data_directory,output_dirctory)

  # Step 2: Call QC plots
  qc_plots <- qc_plots(data_dir_result,output_directory)

  # Step 3: Call filtering function
  filtered_data <- filtering(data_dir_result,output_directory)

  # Step 4: Call sctransform integration
  sctransform_data <- sctransform_integration(filtered_data,output_directory)

  # Step 5: Call clustering
  pca_clusters_result <- pca_clustering(sctransform_data,output_directory)

  # Step 6: Call qc_plots function
  tweedieverse_result <- tweedieverse(pca_clusters_result,output_directory)

  # Combine results into a list
  results_list <- list(
    data_directory = data_dir_result,
    filtering = filtered_data,
    pca_clusters = pca_clusters_result,
    qc_plots = qc_plots_result
  )

  return(results_list)
}
