#' Construct taxon adjacency using phylogenetic distance
#'
#' Computes pairwise cophenetic distances among taxa from the phylogenetic tree
#' in a phyloseq object and converts the distance matrix into a binary
#' taxon-level adjacency matrix.
#'
#'@details
#' This function uses the phylogenetic tree stored in a phyloseq object to define
#' taxon-level dependence. Pairwise distances are computed using cophenetic
#' distances, which represent the branch-length distance between taxa at the tips
#' of the tree. Smaller cophenetic distances indicate that two taxa are more
#' closely related evolutionarily.
#'
#' The resulting distance matrix is converted into a binary adjacency matrix
#' that can be used as the taxon-level neighborhood matrix in SpaMixed. The
#' adjacency matrix can be constructed in three ways:
#'
#' \enumerate{
#'   \item If `threshold` is provided, taxa with cophenetic distance less than or
#'   equal to `threshold` are connected.
#'   \item If `top_k` is provided, each taxon is connected to its top-k closest
#'   taxa in phylogenetic distance before symmetrization.
#'   \item If neither `threshold` nor `top_k` is provided, the threshold is set
#'   to the `quantile_prob` quantile of the off-diagonal pairwise distances. For
#'   example, `quantile_prob = 0.10` connects taxon pairs whose distances are
#'   among the lowest 10% of all pairwise phylogenetic distances.
#' }
#'
#' When `top_k` is used, `sym_rule` controls how the directed k-nearest-neighbor
#' graph is converted to an undirected adjacency matrix. If `sym_rule = "union"`,
#' two taxa are connected if either taxon is among the other's top-k closest
#' taxa. If `sym_rule = "mutual"`, two taxa are connected only if each taxon is
#' among the other's top-k closest taxa. The `"union"` rule is less restrictive
#' and usually produces a denser graph, whereas `"mutual"` is stricter and
#' produces a sparser graph. The `sym_rule` argument is ignored when a
#' threshold-based adjacency matrix is used.
#'
#' If `scale_distance = TRUE`, the cophenetic distances are divided by the
#' median positive pairwise distance. This can be useful when distances are on an
#' arbitrary branch-length scale, but it does not change the relative ordering of
#' distances.
#' 
#' @param phy_obj A phyloseq object containing an OTU table and phylogenetic tree.
#' @param threshold Optional numeric distance threshold. If provided, taxa with distance less than or
#'   equal to this value are connected.
#' @param quantile_prob Numeric distance quantile used to define threshold when `threshold` is NULL and `top_k`
#'   is NULL. Default is 0.10, corresponding to connecting taxon pairs among the closest 10% of pairwise phylogenetic distances.
#' @param top_k Optional integer. If provided, each taxon is connected to its
#'   top-k closest taxa before symmetrization.
#' @param sym_rule Symmetrization rule when `top_k` is provided. Either "union" or "mutual". Default is "union".
#' @param scale_distance Logical. If TRUE, divides distances by their median
#'   positive pairwise distance. Default is FALSE.
#'   
#' @return A list containing:
#' \describe{
#'   \item{distance}{Taxon-by-taxon phylogenetic distance matrix used to construct the adjacency matrix.}
#'   \item{raw_distance}{Original cophenetic distance matrix before optional scaling.}
#'   \item{adj_micro}{Binary taxon adjacency matrix.}
#'   \item{method}{Description of the similarity method.}
#' }
#'
#'
#' @seealso 
#' \code{\link{cooccurrence_similarity}},
#' \code{\link{realized_niche_similarity}},
#' \code{\link{environmental_niche_similarity}},
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
#'phylo_res <- phylogenetic_similarity(
#'phy_obj = physeq_example,
#'quantile_prob = 0.10
#')
#'adj_micro_phylo <- phylo_res$adj_micro
#' dim(adj_micro_phylo)
#' 
#' @export
phylogenetic_similarity <- function(phy_obj,
                                    threshold = NULL,
                                    quantile_prob = 0.10,
                                    top_k = NULL,
                                    sym_rule = c("union", "mutual"),
                                    scale_distance = FALSE) {
  
  sym_rule <- match.arg(sym_rule)
  
  if (!inherits(phy_obj, "phyloseq")) {
    stop("phy_obj must be a phyloseq object.")
  }
  
  if (is.null(phyloseq::phy_tree(phy_obj))) {
    stop("phy_obj must contain a phylogenetic tree.")
  }
  
  tree <- phyloseq::phy_tree(phy_obj)
  taxa_ids <- phyloseq::taxa_names(phy_obj)
  
  if (is.null(tree$edge.length)) {
    warning(
      "The phylogenetic tree has no branch lengths. ",
      "Cophenetic distances may reflect topology rather than branch-length distances."
    )
  }
  
  # --------------------------------------------------
  # Align tree tips with taxa in the phyloseq object
  # --------------------------------------------------
  missing_taxa <- setdiff(taxa_ids, tree$tip.label)
  
  if (length(missing_taxa) > 0) {
    stop(
      "Some taxa in phy_obj are not present in the phylogenetic tree: ",
      paste(missing_taxa, collapse = ", ")
    )
  }
  
  extra_tips <- setdiff(tree$tip.label, taxa_ids)
  
  if (length(extra_tips) > 0) {
    tree <- ape::drop.tip(tree, extra_tips)
  }
  
  # --------------------------------------------------
  # 1. Cophenetic phylogenetic distance matrix
  # --------------------------------------------------
  D_raw <- ape::cophenetic.phylo(tree)
  D_raw <- D_raw[taxa_ids, taxa_ids]
  
  stopifnot(all(rownames(D_raw) == taxa_ids))
  stopifnot(all(colnames(D_raw) == taxa_ids))
  
  D <- D_raw
  
  if (scale_distance) {
    pos_D <- D[upper.tri(D)]
    pos_D <- pos_D[is.finite(pos_D) & pos_D > 0]
    scale_val <- stats::median(pos_D, na.rm = TRUE)
    
    if (!is.na(scale_val) && scale_val > 0) {
      D <- D / scale_val
    }
  }
  
  p <- nrow(D)
  
  # --------------------------------------------------
  # 2. Construct binary adjacency matrix from distance
  # --------------------------------------------------
  R_micro <- matrix(
    0,
    nrow = p,
    ncol = p,
    dimnames = list(taxa_ids, taxa_ids)
  )
  
  if (!is.null(top_k)) {
    
    # k-nearest-neighbor graph in phylogenetic distance space
    for (i in seq_len(p)) {
      d_i <- D[i, ]
      d_i[i] <- Inf
      
      k_i <- min(top_k, p - 1)
      nn_idx <- order(d_i, decreasing = FALSE)[seq_len(k_i)]
      
      R_micro[i, nn_idx] <- 1
    }
    
    # Symmetrize kNN graph
    if (sym_rule == "union") {
      R_micro <- ((R_micro + t(R_micro)) > 0) * 1
    } else {
      R_micro <- ((R_micro + t(R_micro)) == 2) * 1
    }
    
    diag(R_micro) <- 0
    distance_threshold <- NA_real_
    
  } else {
    
    # Quantile-threshold graph
    if (is.null(threshold)) {
      dvec <- D[upper.tri(D)]
      dvec <- dvec[is.finite(dvec)]
      distance_threshold <- stats::quantile(dvec, probs = quantile_prob, na.rm = TRUE)
    } else {
      distance_threshold <- threshold
    }
    
    R_micro <- (D <= distance_threshold) * 1
    diag(R_micro) <- 0
    
    # Ensure symmetry
    R_micro <- pmax(R_micro, t(R_micro))
  }
  
  # --------------------------------------------------
  # 4. Return results
  # --------------------------------------------------
  return(list(
    distance = D,
    raw_distance = D_raw,
    adj_micro = R_micro,
    method = "phylogenetic cophenetic distance"
  ))
}