# library(phyloseq)
# library(spdep)
############################################################
## Region similarity / neighborhood matrix
## Helper 0: Prepare adjacency matrix
############################################################
region_adjacency_matrix <- function(adj, ids = NULL, symmetrize = TRUE) {
  
  if (inherits(adj, "nb")) {
    if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is required when adj is an nb object.")
    }
    library(spdep)
    adj_mat <- spdep::nb2mat(adj, style = "B", zero.policy = TRUE)
    
    if (is.null(ids)) {
      ids <- attr(adj, "region.id")
    }
    
  } else {
    adj_mat <- as.matrix(adj)
  }
  
  storage.mode(adj_mat) <- "numeric"
  
  if (!is.null(ids)) {
    ids <- as.character(ids)
    
    if (nrow(adj_mat) != length(ids)) {
      stop("Length of ids must match number of rows in adj.")
    }
    
    rownames(adj_mat) <- ids
    colnames(adj_mat) <- ids
  }
  
  if (is.null(rownames(adj_mat)) || is.null(colnames(adj_mat))) {
    stop("Adjacency matrix must have row and column names, or ids must be provided.")
  }
  
  diag(adj_mat) <- 0
  
  if (symmetrize) {
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  return(adj_mat)
}

############################################################
## Taxon similarity / neighborhood matrix utilities
## Input: OTU table, preferably taxa x samples
############################################################
############################################################
## Helper 1: convert phyloseq OTU table to taxa x samples matrix
############################################################
as_taxa_by_samples <- function(otu, taxa_are_rows = TRUE) {
  
  # Accept phyloseq object, otu_table, matrix, or data.frame
  if (inherits(otu, "phyloseq")) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("The phyloseq package is required for phyloseq input.")
    } 
    library(phyloseq)
    mat <- as(phyloseq::otu_table(otu), "matrix")
    if (!phyloseq::taxa_are_rows(otu)) {
      mat <- t(mat)
    }
  } else if (inherits(otu, "otu_table")) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("The phyloseq package is required for otu_table input.")
    }
    mat <- as(otu, "matrix")
    if (!phyloseq::taxa_are_rows(otu)) {
      mat <- t(mat)
    }
  } else {
    mat <- as.matrix(otu)
    if (!taxa_are_rows) {
      mat <- t(mat)
    }
  }
  
  mat <- matrix(
    as.numeric(mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  # Remove empty taxa and empty samples
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  
  return(mat)
}


############################################################
## Helper 2: CLR transformation
############################################################
clr_transform <- function(count_mat, pseudocount = 0.5) {
  
  count_mat <- as.matrix(count_mat)
  count_pc <- count_mat + pseudocount
  
  rel_mat <- sweep(count_pc, 2, colSums(count_pc), FUN = "/")
  log_rel <- log(rel_mat)
  
  # CLR within each sample
  clr_mat <- sweep(log_rel, 2, colMeans(log_rel), FUN = "-")
  
  return(clr_mat)  # taxa x samples
}


############################################################
## Helper 3: construct binary adjacency from similarity
############################################################
adjacency_from_similarity <- function(S,
                                      threshold = NULL,
                                      quantile_prob = 0.90,
                                      top_k = NULL,
                                      sym_rule = c("union", "mutual")) {
  
  sym_rule <- match.arg(sym_rule)
  
  S <- as.matrix(S)
  S[is.na(S)] <- 0
  diag(S) <- 0
  
  p <- nrow(S)
  A <- matrix(0, p, p, dimnames = dimnames(S))
  
  if (!is.null(top_k)) {
    # k-nearest-neighbor style adjacency in similarity space
    for (i in seq_len(p)) {
      vals <- S[i, ]
      vals[i] <- -Inf
      
      k_i <- min(top_k, sum(is.finite(vals) & vals > 0))
      if (k_i > 0) {
        idx <- order(vals, decreasing = TRUE)[seq_len(k_i)]
        A[i, idx] <- 1
      }
    }
    
    if (sym_rule == "union") {
      A <- ((A + t(A)) > 0) * 1
    } else {
      A <- ((A + t(A)) == 2) * 1
    }
    
  } else {
    # Threshold-based adjacency
    if (is.null(threshold)) {
      upper_vals <- S[upper.tri(S)]
      upper_vals <- upper_vals[is.finite(upper_vals)]
      threshold <- stats::quantile(upper_vals, probs = quantile_prob, na.rm = TRUE)
    }
    
    A <- (S >= threshold) * 1
    A <- ((A + t(A)) > 0) * 1
  }
  
  diag(A) <- 0
  
  return(A)
}

############################################################
## Helper 4: compute z-scores
############################################################
z_from_marginal <- function(marg) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required.")
  }
  library(INLA)
  mean_val <- inla.emarginal(function(x) x, marg)
  var_val  <- inla.emarginal(function(x) x^2, marg) - mean_val^2
  sd_val   <- sqrt(var_val)
  mean_val / sd_val
}