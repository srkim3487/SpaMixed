#' Construct taxon adjacency using co-occurrence or co-abundance similarity
#'
#' Constructs a taxon-level similarity matrix using CLR-transformed abundance
#' correlations and converts it into a binary adjacency matrix.
#'
#' @details
#' The OTU table is first transformed using a centered log-ratio (CLR)
#' transformation after adding a pseudocount. Pairwise taxon correlations are
#' then computed across samples. If `use_abs_correlation = TRUE`, the absolute
#' value of the correlation is used, so both positive co-occurrence and negative
#' co-exclusion patterns contribute to the similarity matrix. If
#' `use_abs_correlation = FALSE`, only positive correlations are retained.
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
#'   If TRUE, rows are taxa and columns are samples.
#' @param pseudocount Numeric pseudocount added before CLR transformation.
#' @param cor_method Correlation method. Either "spearman" or "pearson".
#' @param use_abs_correlation Logical. If TRUE, uses absolute correlations,
#'   capturing both co-occurrence and co-exclusion strength. If FALSE, only
#'   positive correlations are used.
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
#'   \item{similarity}{Taxon-by-taxon similarity matrix.}
#'   \item{correlation}{Taxon-by-taxon CLR correlation matrix.}
#'   \item{adj_micro}{Binary taxon adjacency matrix.}
#'   \item{counts}{Taxa-by-samples count matrix used in the calculation.}
#'   \item{method}{Description of the similarity method.}
#' }
#'
#' @seealso 
#' \code{\link{realized_niche_similarity}},
#' \code{\link{environmental_niche_similarity}},
#' \code{\link{phylogenetic_similarity}},
#' \code{\link{region_adjacency_matrix}},
#' \code{\link{fit_SpaMixed}}
#' 
#' 
#' @examples
#' 
#' library(SpaMixed)
#' 
#' # Load example data included in the package
#' data(physeq_example)
#'
#' cooc_res <- cooccurrence_similarity(
#'   otu = physeq_example,
#'   taxa_are_rows = TRUE,
#'   pseudocount = 0.5,
#'   cor_method = "spearman",
#'   use_abs_correlation = TRUE,
#'   top_k = 5,
#'   sym_rule = "union"
#' )
#'
#' adj_micro_cooc <- cooc_res$adj_micro
#' dim(adj_micro_cooc)
#' 
#' @export
cooccurrence_similarity <- function(otu,
                                    taxa_are_rows = TRUE,
                                    pseudocount = 0.5,
                                    cor_method = c("spearman", "pearson"),
                                    use_abs_correlation = TRUE,
                                    threshold = NULL,
                                    quantile_prob = 0.90,
                                    top_k = NULL,
                                    sym_rule = c("union", "mutual")) {
  
  cor_method <- match.arg(cor_method)
  sym_rule <- match.arg(sym_rule)
  
  count_mat <- as_taxa_by_samples(otu, taxa_are_rows = taxa_are_rows)
  clr_mat <- clr_transform(count_mat, pseudocount = pseudocount)
  
  # cor() expects variables as columns, so transpose: samples x taxa
  C <- stats::cor(
    t(clr_mat),
    method = cor_method,
    use = "pairwise.complete.obs"
  )
  
  C[is.na(C)] <- 0
  diag(C) <- 0
  
  if (use_abs_correlation) {
    # Captures both co-occurrence and co-exclusion strength
    S <- abs(C)
  } else {
    # Captures only positive co-occurrence
    S <- pmax(C, 0)
  }
  
  A <- adjacency_from_similarity(
    S,
    threshold = threshold,
    quantile_prob = quantile_prob,
    top_k = top_k,
    sym_rule = sym_rule
  )
  
  return(list(
    similarity = S,
    correlation = C,
    adj_micro = A,
    counts = count_mat,
    method = "CLR correlation co-occurrence"
  ))
}
