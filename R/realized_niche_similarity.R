#' Construct taxon adjacency using realized niche overlap
#'
#' Constructs a taxon-level similarity matrix based on how similarly taxa are
#' distributed across samples and converts it into a binary taxon-level adjacency
#' matrix.
#'
#' @details
#' This function treats samples as observed ecological environments and compares
#' taxa according to their abundance distributions across samples. For each taxon,
#' its abundance profile across samples is normalized to sum to one, yielding a
#' sample-distribution profile. Pairwise taxon similarity is then computed using
#' a Pianka-like niche overlap index:
#'
#' \deqn{
#' S_{ik} =
#' \frac{\sum_s P_{is} P_{ks}}
#' {\sqrt{\sum_s P_{is}^2} \sqrt{\sum_s P_{ks}^2}},
#' }
#'
#' where \eqn{P_{is}} is the normalized abundance of taxon \eqn{i} across sample
#' \eqn{s}. Larger values indicate that two taxa tend to occur across similar
#' samples. This measure is an OTU-table-based proxy for realized niche overlap.
#'
#' If `use_relative_abundance = TRUE`, counts are first converted to sample-wise
#' relative abundances before computing taxon distributions across samples. This
#' is often appropriate for microbiome count data because sequencing depths vary
#' across samples. If `use_relative_abundance = FALSE`, raw counts are used.
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
#' taxa. If `sym_rule = "mutual"`, two taxa are connected only if each taxon is
#' among the other's top-k most similar taxa. The `"union"` rule is less
#' restrictive and usually gives a denser graph, whereas `"mutual"` is stricter
#' and gives a sparser graph. The `sym_rule` argument is ignored when a
#' threshold-based adjacency matrix is used.
#'
#' @param otu A phyloseq object, otu_table, matrix, or data.frame.
#' @param taxa_are_rows Logical. Used only when `otu` is a matrix or data.frame.
#' If TRUE, rows are taxa and columns are samples.
#' @param use_relative_abundance Logical. If TRUE, converts counts to relative
#'   abundances before computing overlap.
#' @param threshold Optional numeric threshold for similarity. If provided, taxa with similarity
#'   greater than or equal to this value are connected.
#' @param quantile_prob Numeric quantile used to define the similarity threshold
#'   when `threshold` is NULL and `top_k` is NULL. Default is 0.9, corresponding
#'   to connecting taxon pairs in the top 10% of pairwise similarities.
#' @param top_k Optional integer. If provided, each taxon is connected to its
#'   top-k most similar taxa before symmetrization.
#' @param sym_rule Symmetrization rule used when `top_k` is provided. Either "union" or "mutual". Default is `"union".
#' 
#' @return A list containing:
#' \describe{
#'   \item{similarity}{Taxon-by-taxon realized niche-overlap similarity matrix.}
#'   \item{adj_micro}{Binary taxon adjacency matrix.}
#'   \item{counts}{Taxa-by-samples count matrix used in the calculation.}
#'   \item{method}{Description of the similarity method.}
#' }
#' 
#'
#' @seealso 
#' \code{\link{cooccurrence_similarity}},
#' \code{\link{environmental_niche_similarity}},
#' \code{\link{phylogenetic_similarity}},
#' \code{\link{region_adjacency_matrix}},
#' \code{\link{fit_SpaMixed}}
#' 
#' @examples
#' 
#' library(SpaMixed)
#' 
#' # Load example data included in the package
#' data(physeq_example)
#'
#'niche_res <- realized_niche_similarity(
#'otu = physeq_example, 
#'taxa_are_rows = TRUE, 
#'top_k = 5)
#'
#' adj_micro_niche <- niche_res$adj_micro
#' dim(adj_micro_niche)
#'
#' @export
realized_niche_similarity <- function(otu,
                                      taxa_are_rows = TRUE,
                                      use_relative_abundance = TRUE,
                                      threshold = NULL,
                                      quantile_prob = 0.90,
                                      top_k = NULL,
                                      sym_rule = c("union", "mutual")) {
  
  count_mat <- as_taxa_by_samples(otu, taxa_are_rows = taxa_are_rows)
  
  if (use_relative_abundance) {
    # Sample-wise relative abundance
    X <- sweep(count_mat, 2, colSums(count_mat), FUN = "/")
  } else {
    X <- count_mat
  }
  
  # For each taxon, convert its abundance across samples into a distribution
  row_totals <- rowSums(X, na.rm = TRUE)
  P <- sweep(X, 1, row_totals, FUN = "/")
  P[is.na(P)] <- 0
  
  # Pianka-like niche overlap:
  # S_ik = sum_s P_is P_ks / sqrt(sum_s P_is^2 * sum_s P_ks^2)
  denom <- sqrt(rowSums(P^2, na.rm = TRUE))
  S <- P %*% t(P)
  S <- S / outer(denom, denom)
  S[is.na(S)] <- 0
  diag(S) <- 1
  
  A <- adjacency_from_similarity(
    S,
    threshold = threshold,
    quantile_prob = quantile_prob,
    top_k = top_k,
    sym_rule = sym_rule
  )
  
  return(list(
    similarity = S,
    adj_micro = A,
    counts = count_mat,
    method = "OTU-only realized niche overlap"
  ))
}
