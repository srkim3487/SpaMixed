#' Prepare a region-level adjacency matrix
#'
#' Converts a region-level neighborhood object into a numeric adjacency matrix
#' for use in spatial random effects.
#'
#' @details
#' This function prepares the region-level neighborhood matrix used to model
#' spatial dependence among regions or sampling locations in SpaMixed. The input
#' can be either an \code{spdep} \code{nb} object or a square adjacency matrix.
#'
#' If `adj` is an \code{nb} object, the function converts it to a binary
#' adjacency matrix using \code{\link[spdep:nb2mat]{spdep::nb2mat}} with
#' `style = "B"`. If `ids` is not provided, the function attempts to use
#' `attr(adj, "region.id")` as the region identifiers.
#'
#' If `adj` is already a matrix, it is converted to a numeric matrix. Row and
#' column names must be present, unless `ids` is provided. If `ids` is provided,
#' it is used to assign both row and column names. The diagonal is set to zero
#' because regions are not treated as neighbors of themselves.
#'
#' When `symmetrize = TRUE`, the adjacency matrix is made symmetric using
#' \code{pmax(adj_mat, t(adj_mat))}. This is useful when the input neighborhood
#' object is directed or when one region lists another as a neighbor but not
#' vice versa. The resulting matrix can be used as the region-level graph in
#' spatial random-effect models, such as the \code{besagproper2} model in
#' \pkg{INLA}.
#'
#' @param adj A region-level neighborhood object. This can be either an
#'   \code{spdep} \code{nb} object or a square matrix-like object representing
#'   region adjacency.
#' @param ids Optional character vector of region identifiers. If provided,
#'   `ids` is used as both row and column names of the adjacency matrix. The
#'   length of `ids` must match the number of regions.
#' @param symmetrize Logical. If TRUE, the adjacency matrix is symmetrized using
#'   \code{pmax(adj_mat, t(adj_mat))}. Default is TRUE.
#'
#' @return A numeric square region-by-region adjacency matrix with zero diagonal
#'   and row and column names corresponding to region identifiers.
#'
#'
#' @seealso
#' \code{\link{cooccurrence_similarity}},
#' \code{\link{realized_niche_similarity}},
#' \code{\link{environmental_niche_similarity}},
#' \code{\link{phylogenetic_similarity}},
#' \code{\link{fit_SpaMixed}}
#' 
#' @examples
#' adj <- matrix(
#'   c(
#'     0, 1, 0,
#'     1, 0, 1,
#'     0, 1, 0
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' region_adj <- region_adjacency_matrix(
#'   adj = adj,
#'   ids = c("region1", "region2", "region3")
#' )
#'
#' region_adj
#' 
#'
#' @export
region_adjacency_matrix <- function(adj, ids = NULL, symmetrize = TRUE) {
  
  if (inherits(adj, "nb")) {
    
    if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is required when `adj` is an nb object.")
    }
    
    adj_mat <- spdep::nb2mat(
      adj,
      style = "B",
      zero.policy = TRUE
    )
    
    if (is.null(ids)) {
      ids <- attr(adj, "region.id")
    }
    
  } else {
    
    adj_mat <- as.matrix(adj)
  }
  
  if (nrow(adj_mat) != ncol(adj_mat)) {
    stop("`adj` must be a square adjacency matrix or an nb object.")
  }
  
  storage.mode(adj_mat) <- "numeric"
  
  if (!is.null(ids)) {
    
    ids <- as.character(ids)
    
    if (nrow(adj_mat) != length(ids)) {
      stop("Length of `ids` must match the number of rows in `adj`.")
    }
    
    rownames(adj_mat) <- ids
    colnames(adj_mat) <- ids
  }
  
  if (is.null(rownames(adj_mat)) || is.null(colnames(adj_mat))) {
    stop(
      "Adjacency matrix must have row and column names, ",
      "or `ids` must be provided."
    )
  }
  
  if (!identical(rownames(adj_mat), colnames(adj_mat))) {
    stop("Row names and column names of the adjacency matrix must match.")
  }
  
  diag(adj_mat) <- 0
  
  if (symmetrize) {
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  return(adj_mat)
}