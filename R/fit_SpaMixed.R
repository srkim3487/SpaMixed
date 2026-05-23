#' Fit the SpaMixed model
#'
#' Fits the spatial mixed model (SpaMixed) for high-dimensional
#' microbiome count data. The model evaluates associations between an exposure
#' variable and microbial taxa while accounting for region-level spatial
#' dependence, sample-level heterogeneity, and taxon-level dependence.
#'
#' @details
#' `fit_SpaMixed()` fits a zero-inflated Poisson mixed model to microbiome count
#' data stored in a \code{\link[phyloseq:phyloseq-class]{phyloseq}} object.
#' The exposure of interest is modeled through taxon-specific exposure effects,
#' allowing each taxon to have its own association with the exposure. Additional
#' covariates can also be adjusted using taxon-specific fixed effects.
#'
#' The model includes three structured components:
#' \enumerate{
#'   \item a region-level spatial random effect, specified using a region
#'   adjacency matrix;
#'   \item a sample-level random effect, which captures sample-specific
#'   heterogeneity;
#'   \item a taxon-level random effect, specified using a taxon adjacency matrix
#'   such as a phylogenetic, co-occurrence, niche-overlap, or other biologically
#'   motivated similarity structure.
#' }
#'
#' Model fitting is performed using \pkg{INLA}. Feature selection is performed
#' by computing posterior Wald statistics for the taxon-specific exposure
#' effects and applying local false discovery rate (lfdr) control using
#' \pkg{fdrtool}. Selected taxa are returned together with posterior summaries
#' of the exposure effects. If a taxonomy table is available in `phy_obj`, the
#' selected features are also annotated with taxonomic information.
#'
#' Continuous variables listed in `scale_vars` are centered and scaled to unit
#' variance before model fitting. Variables listed in `factor_vars` are converted
#' to factors before constructing the model matrix. The exposure variable must be
#' numeric; categorical exposures should be converted to dummy variables before
#' calling this function.
#'
#' @param phy_obj A \code{\link[phyloseq:phyloseq-class]{phyloseq}} object
#'   containing an OTU/count table, sample metadata, and optionally a taxonomy
#'   table and phylogenetic tree.
#' @param exposure Character string giving the name of the exposure variable in
#'   \code{sample_data(phy_obj)}. This variable is modeled as the primary
#'   exposure of interest and must be numeric.
#' @param covariates Optional character vector giving the names of additional
#'   covariates in \code{sample_data(phy_obj)} to adjust for. The exposure should
#'   not be included in `covariates`.
#' @param region_var Character string giving the name of the region or spatial
#'   unit variable in \code{sample_data(phy_obj)}. This variable is used to link
#'   samples to the region-level adjacency matrix.
#' @param region_adj A binary region-by-region adjacency matrix describing
#'   spatial neighborhood relationships among regions. Row and column names
#'   should correspond to region identifiers.
#' @param region_ids Optional character vector of region identifiers. This
#'   argument is currently retained for compatibility and can be used to document
#'   the intended ordering of regions.
#' @param micro_adj A binary taxon-by-taxon adjacency matrix describing
#'   neighborhood relationships among taxa. Row and column names should match
#'   the taxon order used by the model. This matrix can be constructed using
#'   functions such as \code{\link{phylogenetic_similarity}},
#'   \code{\link{cooccurrence_similarity}},
#'   \code{\link{realized_niche_similarity}}, or
#'   \code{\link{environmental_niche_similarity}}.
#' @param scale_vars Optional character vector of continuous variables in
#'   \code{sample_data(phy_obj)} to center and scale before model fitting.
#' @param factor_vars Optional character vector of variables in
#'   \code{sample_data(phy_obj)} to convert to factors before model fitting.
#' @param num_threads Integer specifying the number of threads passed to
#'   \code{\link[INLA:inla]{INLA::inla}}. Default is 2.
#' @param control_fixed A list passed to the `control.fixed` argument of
#'   \code{\link[INLA:inla]{INLA::inla}}. Default is
#'   \code{list(prec.intercept = 0.001)}.
#' @param control_compute A list passed to the `control.compute` argument of
#'   \code{\link[INLA:inla]{INLA::inla}}. Default is
#'   \code{list(config = TRUE, dic = TRUE, waic = TRUE)}.
#' @param exposure_prec_prior Numeric vector of length 2 specifying the
#'   parameters of the PC prior for the precision of the taxon-specific exposure
#'   effects. Default is \code{c(0.5, 0.01)}.
#' @param region_prec_prior Numeric vector of length 2 specifying the PC prior
#'   parameters for the precision of the region-level spatial random effect.
#'   Default is \code{c(2, 0.05)}.
#' @param region_lambda_prior Numeric vector of length 2 specifying the lbeta
#'   prior parameters for the spatial dependence parameter in the
#'   region-level \code{besagproper2} model. Default is \code{c(1, 1)}.
#' @param sample_prec_prior Numeric vector of length 2 specifying the PC prior
#'   parameters for the precision of the sample-level random effect. Default is
#'   \code{c(2, 0.05)}.
#' @param micro_prec_prior Numeric vector of length 2 specifying the PC prior
#'   parameters for the precision of the taxon-level random effect. Default is
#'   \code{c(2, 0.05)}.
#' @param micro_lambda_prior Numeric vector of length 2 specifying the beta
#'   prior parameters for the dependence parameter in the taxon-level
#'   \code{besagproper2} model. Default is \code{c(2, 2)}.
#' @param lfdr_thres Numeric local false discovery rate threshold used for
#'   feature selection. Default is 0.2.
#'
#' @return An object of class `"SpaMixedFit"`, which is a list containing:
#' \describe{
#'   \item{fit}{The fitted \pkg{INLA} model object.}
#'   \item{fixed_effect}{Posterior summaries of fixed effects from the fitted model.}
#'   \item{exposure_effect}{Posterior summaries of taxon-specific exposure effects.}
#'   \item{random_effect}{Posterior summaries of model hyperparameters.}
#'   \item{selected_feature}{A data frame of taxa selected by the lfdr criterion,
#'   including posterior summaries and, when available, taxonomy annotations.}
#' }
#'
#' @seealso
#' \code{\link{cooccurrence_similarity}},
#' \code{\link{realized_niche_similarity}},
#' \code{\link{environmental_niche_similarity}},
#' \code{\link{phylogenetic_similarity}},
#' \code{\link{region_adjacency_matrix}}
#' 
#' @examples
#' library(SpaMixed)
#' library(spdep)
#'
#' data(physeq_example)
#'
#' # Construct a taxon-level adjacency matrix using co-occurrence similarity
#' cooc_res <- cooccurrence_similarity(
#' otu = physeq_example,
#' taxa_are_rows = TRUE,
#' pseudocount = 0.5,
#' cor_method = "spearman",
#' use_abs_correlation = TRUE,
#' top_k = 5,
#' sym_rule = "union"
#' )
#' adj_micro <- cooc_res$adj_micro
#'
#' # Construct a region-level adjacency matrix 
#' # As an example, we used a 4-nearest neighborhood structure with a torus on a 9*9 regular lattice.
#' # region_adj should be a binary region-by-region adjacency matrix
#' # with row and column names matching the region variable in sample_data.
#' 
#' nbd_index <- cell2nb(nrow = 9, ncol = 9, type = "rook", torus = TRUE) 
#' adj_mat <- region_adjacency_matrix(nbd_index)
#' attr(adj_mat, "region.id") <- 1:81
#' rownames(adj_mat) <- colnames(adj_mat) <- attr(adj_mat, "region.id") 
#' 
#' # Fit SpaMixed 
#' fit <- fit_SpaMixed(
#'   phy_obj = physeq_example,
#'   exposure = "PM2.5",
#'   covariates = c("age", "gender"),
#'   region_var = "zip_code",
#'   region_adj = adj_mat,
#'   micro_adj = adj_micro,
#'   scale_vars = c("PM2.5", "age"),
#'   factor_vars = c("gender"),
#'   lfdr_thres = 0.2
#' )
#'
#' fit$selected_feature
#' 
#'
#' @export
fit_SpaMixed <- function(
    phy_obj,                   
    exposure,
    covariates = NULL,
    region_var,
    region_adj,
    region_ids = NULL,
    micro_adj,

    scale_vars = NULL,
    factor_vars = NULL,
    
    num_threads = 10,
    control_fixed = list(prec.intercept = 0.001),
    control_compute = list(config = TRUE, dic = TRUE, waic = TRUE),
    exposure_prec_prior = c(0.5, 0.01),
    region_prec_prior = c(2, 0.05),
    region_lambda_prior = c(1, 1),
    sample_prec_prior = c(2, 0.05),
    micro_prec_prior = c(2, 0.05),
    micro_lambda_prior = c(2, 2),

    lfdr_thres = 0.2
) {
  
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required.")
  }
  
  if (!inherits(phy_obj, "phyloseq")) {
    stop("phy_obj must be a phyloseq object.")
  }
  
  library(phyloseq)
  library(ape)
  library(INLA)

  
  if (!is.null(covariates) && exposure %in% covariates) {
    stop("Do not include exposure in covariates. The exposure is modeled separately.")
  }
  
  ##########################################################
  ## 1. Extract counts and sample metadata
  ##########################################################
  count_mat <- as_taxa_by_samples(phy_obj)
  
  taxa_ids <- rownames(count_mat)
  sample_ids <- colnames(count_mat)
  
  J <- nrow(count_mat)
  N <- ncol(count_mat)
  
  sample_dt <- as(phyloseq::sample_data(phy_obj), "data.frame") 
  
  if (!all(sample_ids %in% rownames(sample_dt))) {
    stop("Sample names in otu_table are not all present in sample_data.")
  }
  
  sample_dt <- sample_dt[sample_ids, , drop = FALSE]
  
  required_vars <- unique(c(exposure, covariates, region_var))
  missing_vars <- setdiff(required_vars, colnames(sample_dt))
  
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from sample_data: ",
         paste(missing_vars, collapse = ", "))
  }
  
  ##########################################################
  ## 2. Make safe metadata names
  ##########################################################
  
  original_names <- colnames(sample_dt)
  safe_names <- make.names(original_names, unique = TRUE)
  name_map <- setNames(safe_names, original_names)
  
  colnames(sample_dt) <- safe_names
  
  exposure_safe <- name_map[[exposure]]
  region_safe <- name_map[[region_var]]
  covariates_safe <- if (!is.null(covariates)) name_map[covariates] else NULL
  scale_vars_safe <- if (!is.null(scale_vars)) name_map[scale_vars] else NULL
  factor_vars_safe <- if (!is.null(factor_vars)) name_map[factor_vars] else NULL
  
  ##########################################################
  ## 3. Scale continuous variables and convert factors
  ##########################################################
  
  if (!is.null(scale_vars_safe)) {
    for (v in scale_vars_safe) {
      sample_dt[[v]] <- as.numeric(scale(sample_dt[[v]], center = TRUE, scale = TRUE))
    }
  }
  
  if (!is.null(factor_vars_safe)) {
    for (v in factor_vars_safe) {
      sample_dt[[v]] <- as.factor(sample_dt[[v]])
    }
  }
  
  sample_dt[[region_safe]] <- as.factor(sample_dt[[region_safe]])
  
  if (!is.numeric(sample_dt[[exposure_safe]])) {
    stop("The exposure variable must be numeric. If it is categorical, create dummy variables first.")
  }
  
  ##########################################################
  ## 4. Remove samples with missing model variables
  ##########################################################
  
  model_vars <- unique(c(exposure_safe, covariates_safe, region_safe))
  keep_samples <- stats::complete.cases(sample_dt[, model_vars, drop = FALSE])
  
  if (!all(keep_samples)) {
    message("Removing ", sum(!keep_samples), " samples with missing model variables.")
    sample_dt <- sample_dt[keep_samples, , drop = FALSE]
    count_mat <- count_mat[, rownames(sample_dt), drop = FALSE]
    sample_ids <- colnames(count_mat)
    N <- ncol(count_mat)
  }
  
  ##########################################################
  ## 5. Prepare region adjacency matrix
  ##########################################################
  
  regional_values <- rownames(region_adj)
  
  observed_regions <- as.character(sample_dt[[region_safe]])
  missing_regions <- setdiff(unique(observed_regions), regional_values)
  
  if (length(missing_regions) > 0) {
    stop("Some sample regions are not found in region_adj: ",
         paste(missing_regions, collapse = ", "))
  }
  
  ##########################################################
  ## 6. Create long-format response and indexing variables
  ##########################################################
  micro_vec <- as.numeric(count_mat)
  
  taxon_code <- rownames(micro_adj)
  data_long <- data.frame(
    micro_vec = micro_vec,
    taxon_int = factor(rep(taxon_code, N), levels = taxon_code),
    beta_idx = factor(rep(taxon_code, N), levels = taxon_code),
    exposure_inla = rep(sample_dt[[exposure_safe]], each = J),
    region = factor(rep(as.character(sample_dt[[region_safe]]), each = J),
                    levels = regional_values),
    sampling_rnd = factor(rep(sample_ids, each = J), levels = sample_ids),
    micro_rnd = factor(rep(taxon_code, N), levels = taxon_code)
  )
  
  ##########################################################
  ## 7. Build taxon-specific fixed effects for covariates
  ##########################################################
  
  covariate_design_cols <- character(0)
  covariate_design_map <- NULL
  
  if (!is.null(covariates_safe) && length(covariates_safe) > 0) {
    
    cov_formula <- stats::as.formula(
      paste("~", paste(covariates_safe, collapse = " + "))
    )
    
    cov_mm <- stats::model.matrix(cov_formula, data = sample_dt)
    
    if ("(Intercept)" %in% colnames(cov_mm)) {
      cov_mm <- cov_mm[, colnames(cov_mm) != "(Intercept)", drop = FALSE]
    }
    
    if (ncol(cov_mm) > 0) {
      
      base_taxon <- stats::model.matrix(
        ~ 0 + taxon_int,
        data = data.frame(taxon_int = data_long$taxon_int)
      )
      
      colnames(base_taxon) <- taxon_code
      
      cov_design_list <- list()
      covariate_design_map <- list()
      
      for (k in seq_len(ncol(cov_mm))) {
        
        x_k <- as.numeric(cov_mm[, k])
        x_long <- rep(x_k, each = J)
        
        X_k <- base_taxon * x_long
        
        cov_name <- make.names(colnames(cov_mm)[k])
        new_colnames <- paste0(cov_name, "_", taxon_code)
        
        colnames(X_k) <- new_colnames
        
        cov_design_list[[k]] <- X_k
        
        covariate_design_map[[k]] <- data.frame(
          original_covariate_column = colnames(cov_mm)[k],
          safe_covariate_column = cov_name,
          taxon_code = taxon_code,
          taxon_id = taxa_ids,
          design_column = new_colnames,
          stringsAsFactors = FALSE
        )
      }
      
      cov_design_mat <- do.call(cbind, cov_design_list)
      covariate_design_cols <- colnames(cov_design_mat)
      
      data_long <- cbind(data_long, as.data.frame(cov_design_mat))
      covariate_design_map <- do.call(rbind, covariate_design_map)
    }
  }
  
  ##########################################################
  ## 8. Build INLA formula
  ##########################################################
  
  formula_terms <- c(
    "0",
    "taxon_int"
  )
  
  if (length(covariate_design_cols) > 0) {
    formula_terms <- c(formula_terms, covariate_design_cols)
  }
  
  
  ## Taxon-specific exposure effect with shrinkage
  exposure_term <- paste0(
    "f(beta_idx, exposure_inla, model = 'iid', ",
    "hyper = list(prec = list(prior = 'pc.prec', param = c(",
    exposure_prec_prior[1], ", ", exposure_prec_prior[2], "))))"
  )
  
  formula_terms <- c(formula_terms, exposure_term)
  
  ## Region-level random effect
  region_term <- paste0(
      "f(region, values = regional_values, model = 'besagproper2', graph = region_adj, ",
      "hyper = list(",
      "prec = list(prior = 'pc.prec', param = c(",
      region_prec_prior[1], ", ", region_prec_prior[2], ")), ",
      "lambda = list(prior = 'logitbeta', param = c(",
      region_lambda_prior[1], ", ", region_lambda_prior[2], "))))"
    )
    
  formula_terms <- c(formula_terms, region_term)
  
  ## Sample-level random effect
  sample_term <- paste0(
      "f(sampling_rnd, model = 'iid', ",
      "hyper = list(prec = list(prior = 'pc.prec', param = c(",
      sample_prec_prior[1], ", ", sample_prec_prior[2], "))))"
    )
    
  formula_terms <- c(formula_terms, sample_term)

  
  ## Taxon-level random effect
  micro_term <- paste0(
      "f(micro_rnd, values = micro_values, model = 'besagproper2', graph = micro_adj, ",
      "hyper = list(",
      "prec = list(prior = 'pc.prec', param = c(",
      micro_prec_prior[1], ", ", micro_prec_prior[2], ")), ",
      "lambda = list(prior = 'logitbeta', param = c(",
      micro_lambda_prior[1], ", ", micro_lambda_prior[2], "))))"
    )
    
  formula_terms <- c(formula_terms, micro_term)
  
  ## Formula
  formula_string <- paste(
    "micro_vec ~",
    paste(formula_terms, collapse = " + ")
  )
  
  formula <- stats::as.formula(formula_string)
  environment(formula) <- environment()
  
  ##########################################################
  ## 9. Fit INLA model
  ##########################################################
  micro_values <- taxon_code

    fit <- INLA::inla(
      formula,
      family = "zeroinflatedpoisson1",
      data = data_long,
      num.threads = num_threads,
      control.fixed = control_fixed,
      control.compute = control_compute
    )
    
    
    ##########################################################
    ## Variable selection
    ##########################################################
    pm_margs_all <- fit$marginals.random$beta_idx
    z_scores <- sapply(pm_margs_all, z_from_marginal)
    
    if (!requireNamespace("fdrtool", quietly = TRUE)) {
      stop("Package 'fdrtool' is required.")
    }
    library(fdrtool)
    
    lfdr_res <- fdrtool(z_scores, statistic = "normal", plot = FALSE, cutoff.method="locfdr", verbose = FALSE)
    sel_stat <- lfdr_res$lfdr      
    sel_idx  <- which(sel_stat < lfdr_thres)
    
    PM_spa <- fit$summary.random$beta_idx
    PM_coef <- cbind(PM_spa, Taxa = rownames(otu_table(phy_obj)))
    plot_pm_taxa <- PM_coef[sel_idx,]
    
    if (!is.null(tax_table(phy_obj, errorIfNULL = FALSE))) {
      tax_df <- as.data.frame(tax_table(phy_obj)) %>%
      tibble::rownames_to_column(var = "taxa_id")
    
      plot_pm_taxa = left_join(plot_pm_taxa, tax_df, by=c("Taxa" = "taxa_id"))
    }
  
  
  
  ##########################################################
  ## 10. Return package-style output
  ##########################################################
  
  out <- list(
    fit = fit,
    fixed_effect = fit$summary.fixed, 
    exposure_effect = fit$summary.random$beta_idx,
    random_effect = fit$summary.hyperpar,
    selected_feature = plot_pm_taxa
  )
  
  class(out) <- "SpaMixedFit"
  
  return(out)
}
