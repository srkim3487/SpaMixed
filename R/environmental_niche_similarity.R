#' Construct taxon adjacency using environmental niche similarity
#'
#' Constructs taxon-level environmental niche centroids using sample-level
#' environmental variables and taxon relative-abundance weights, then computes
#' taxon similarity based on distances between these centroids.
#'
#' @details
#' This function defines taxon-level similarity by comparing the environmental
#' conditions of the samples in which each taxon is observed. First, the OTU
#' table is converted to sample-wise relative abundance. For each taxon, an
#' environmental centroid is computed as a weighted average of the sample-level
#' environmental variables, where the weights are the relative abundances of
#' that taxon across samples. Conceptually, taxa with similar environmental
#' centroids are treated as occupying similar observed, or realized,
#' environmental niches.
#'
#' Specifically, let \eqn{r_{is}} denote the relative abundance of taxon
#' \eqn{i} in sample \eqn{s}, and let \eqn{E_s} denote the environmental
#' covariate vector for sample \eqn{s}. The environmental centroid for taxon
#' \eqn{i} is computed as
#'
#' \deqn{
#' c_i =
#' \frac{\sum_s r_{is} E_s}{\sum_s r_{is}}.
#' }
#'
#' Pairwise Euclidean distances are then computed between taxon centroids.
#' These distances are converted into similarities using a radial basis function
#' (RBF) kernel:
#'
#' \deqn{
#' S_{ik} =
#' \exp\left\{-\frac{D_{ik}^2}{2h^2}\right\},
#' }
#'
#' where \eqn{D_{ik}} is the Euclidean distance between the environmental
#' centroids of taxa \eqn{i} and \eqn{k}, and \eqn{h} is the bandwidth. If
#' `bandwidth` is not provided, the median positive pairwise centroid distance
#' is used.
#'
#' This measure should be interpreted as an observed environmental niche
#' similarity based on the available samples and environmental variables. It is
#' not a full ecological niche model and does not imply causal environmental
#' preference.
#'
#' The binary adjacency matrix can be constructed in three ways:
#' \enumerate{
#'   \item If `threshold` is provided, taxa with similarity greater than or equal
#'   to `threshold` are connected.
#'   \item If `top_k` is provided, each taxon is connected to its top-k most
#'   similar taxa before symmetrization.
#'   \item If neither `threshold` nor `top_k` is provided, the threshold is set
#'   to the `quantile_prob` quantile of the off-diagonal similarity values. For
#'   example, `quantile_prob = 0.90` connects taxa whose similarity is in the
#'   top 10% of all pairwise similarities.
#' }
#'
#' When `top_k` is used, `sym_rule` controls how the directed top-k graph is
#' converted to an undirected adjacency matrix. If `sym_rule = "union"`, two
#' taxa are connected if either taxon is among the other's top-k most similar
#' taxa. If `sym_rule = "mutual"`, two taxa are connected only if both taxa are
#' among each other's top-k most similar taxa. The `"union"` rule is less
#' restrictive and usually gives a denser graph, whereas `"mutual"` is stricter
#' and gives a sparser graph. The `sym_rule` argument is ignored when a
#' threshold-based adjacency matrix is used.
#'
#' Environmental variables can be numeric or categorical. Categorical variables
#' are automatically expanded using `model.matrix()`. If `scale_env = TRUE`,
#' the resulting environmental design matrix is standardized before computing
#' taxon centroids.
#'
#' @param otu A phyloseq object, otu_table, matrix, or data.frame.
#' @param env_df Data.frame of sample-level environmental variables. Row names
#'   must be sample IDs matching the OTU table column names.
#' @param env_cols Character vector of environmental variables to use. If NULL,
#'   all columns in `env_df` are used.
#' @param taxa_are_rows Logical. Used only when `otu` is a matrix or data.frame.
#'   If TRUE, rows are taxa and columns are samples.
#' @param scale_env Logical. If TRUE, environmental variables are standardized
#'   after numeric variables and dummy-coded categorical variables are assembled
#'   into a design matrix. Default is TRUE.
#' @param bandwidth Optional numeric bandwidth for the RBF similarity kernel. If
#'   NULL, the median positive pairwise centroid distance is used.
#' @param threshold Optional numeric similarity threshold. If provided, taxa with
#'   similarity greater than or equal to this value are connected.
#' @param quantile_prob Numeric quantile used to define the similarity threshold
#'   when `threshold` is NULL and `top_k` is NULL. Default is 0.90, corresponding
#'   to connecting taxon pairs in the top 10% of pairwise similarities.
#' @param top_k Optional integer. If provided, each taxon is connected to its
#'   top-k most similar taxa before symmetrization.
#' @param sym_rule Symmetrization rule used when `top_k` is provided. Either
#'   `"union"` or `"mutual"`. Default is `"union"`. 
#'
#' @return A list containing:
#' \describe{
#'   \item{similarity}{Taxon-by-taxon environmental niche similarity matrix.}
#'   \item{distance}{Taxon-by-taxon Euclidean distance matrix between environmental centroids.}
#'   \item{adj_micro}{Binary taxon adjacency matrix.}
#'   \item{centroid}{Taxon-by-environmental-feature centroid matrix.}
#'   \item{counts}{Taxa-by-samples count matrix used in the calculation.}
#'   \item{method}{Description of the similarity method.}
#' }
#'
#'
#' @seealso 
#' \code{\link{cooccurrence_similarity}},
#' \code{\link{realized_niche_similarity}},
#' \code{\link{phylogenetic_similarity}},
#' \code{\link{region_adjacency_matrix}},
#' \code{\link{fit_SpaMixed}}
#' 
#' @examples
#' 
#' library(SpaMixed)
#' library(phyloseq)
#' 
#' # Load example data included in the package
#' data(physeq_example)
#'
#' env_df <- data.frame(PM25=sample_data(physeq_example)$`PM2.5`)
#' rownames(env_df) <- colnames(otu_table(physeq_example))
#' env_res <- environmental_niche_similarity(
#'  otu = physeq_example,
#'  env_df = env_df,
#' env_cols = c("PM25"),
#'  taxa_are_rows = TRUE,
#'  top_k = 5,
#'  sym_rule = "union"
#'  )
#'
#' adj_micro_env_niche <- env_res$adj_micro
#' dim(adj_micro_env_niche)
#' 
#' @export
environmental_niche_similarity <- function(otu,
                                           env_df,
                                           env_cols = NULL,
                                           taxa_are_rows = TRUE,
                                           scale_env = TRUE,
                                           bandwidth = NULL,
                                           threshold = NULL,
                                           quantile_prob = 0.90,
                                           top_k = NULL,
                                           sym_rule = c("union", "mutual")) {
  
  count_mat <- as_taxa_by_samples(otu, taxa_are_rows = taxa_are_rows)
  
  if (is.null(rownames(env_df))) {
    stop("env_df must have sample IDs as row names.")
  }
  
  common_samples <- intersect(colnames(count_mat), rownames(env_df))
  if (length(common_samples) < 3) {
    stop("Too few overlapping samples between OTU table and env_df.")
  }
  
  count_mat <- count_mat[, common_samples, drop = FALSE]
  env_df <- env_df[common_samples, , drop = FALSE]
  
  if (is.null(env_cols)) {
    env_cols <- colnames(env_df)
  }
  
  E <- as.matrix(env_df[, env_cols, drop = FALSE])
  E <- matrix(as.numeric(E), nrow = nrow(E), dimnames = dimnames(E))
  
  keep_samples <- rowSums(is.na(E)) == 0
  E <- E[keep_samples, , drop = FALSE]
  count_mat <- count_mat[, rownames(E), drop = FALSE]
  
  if (scale_env) {
    E <- scale(E)
  }
  
  # Use relative abundance as weights
  rel_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/")
  taxon_weight <- rowSums(rel_mat, na.rm = TRUE)
  
  # Taxon environmental centroid
  centroid <- rel_mat %*% E
  centroid <- sweep(centroid, 1, taxon_weight, FUN = "/")
  centroid[is.na(centroid)] <- 0
  
  D <- as.matrix(stats::dist(centroid, method = "euclidean"))
  
  if (is.null(bandwidth)) {
    pos_D <- D[upper.tri(D)]
    pos_D <- pos_D[is.finite(pos_D) & pos_D > 0]
    bandwidth <- stats::median(pos_D, na.rm = TRUE)
  }
  
  # RBF similarity from environmental niche distance
  S <- exp(-(D^2) / (2 * bandwidth^2))
  diag(S) <- 1
  
  A <- adjacency_from_similarity(
    S,
    threshold = threshold,
    quantile_prob = quantile_prob,
    top_k = top_k
  )
  
  return(list(
    similarity = S,
    distance = D,
    adj_micro = A,
    centroid = centroid,
    counts = count_mat,
    method = "environmental niche similarity"
  ))
}