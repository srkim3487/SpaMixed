This R package implements the SpaMixed method introduced in:

Kim, S., Wang, C., Kwak, S., Bain, A., Segal, L., Ahn, J., and Li, H. Spatial mixed models for assessing environmental exposure effects on the microbiome. 

To install and load the package:
```r
devtools::install_github("srkim3487/SpaMixed")
library(SpaMixed)
```

Below are examples using the included simulated dataset. For details on each function, please refer to the corresponding help page.
```r
library(SpaMixed)
library(spdep)

data(physeq_example)

# Construct taxon adjacency using co-occurrence or co-abundance similarity
cooc_res <- cooccurrence_similarity(
   otu = physeq_example,
   taxa_are_rows = TRUE,
   pseudocount = 0.5,
   cor_method = "spearman",
   use_abs_correlation = TRUE,
   top_k = 5,
   sym_rule = "union"
 )
adj_micro <- cooc_res$adj_micro

# Construct taxon adjacency using realized niche overlap
niche_res <- realized_niche_similarity(
otu = physeq_example, 
taxa_are_rows = TRUE, 
top_k = 5)

adj_micro_niche <- niche_res$adj_micro

# Construct taxon adjacency using environmental niche similarity
env_df <- data.frame(PM25=sample_data(physeq_example)$`PM2.5`)
rownames(env_df) <- colnames(otu_table(physeq_example))
env_res <- environmental_niche_similarity(
  otu = physeq_example,
  env_df = env_df,
 env_cols = c("PM25"),
  taxa_are_rows = TRUE,
  top_k = 5,
  sym_rule = "union"
  )

adj_micro_env_niche <- env_res$adj_micro

# Construct taxon adjacency using phylogenetic distance
phylo_res <- phylogenetic_similarity(
phy_obj = physeq_example,
quantile_prob = 0.10
)
adj_micro_phylo <- phylo_res$adj_micro

# Construct a region-level adjacency matrix 
# As an example, we used a 4-nearest neighborhood structure with a torus on a 9*9 regular lattice.
# region_adj should be a binary region-by-region adjacency matrix
# with row and column names matching the region variable in sample_data.

nbd_index <- cell2nb(nrow = 9, ncol = 9, type = "rook", torus = TRUE) 
adj_mat <- region_adjacency_matrix(nbd_index)
attr(adj_mat, "region.id") <- 1:81
rownames(adj_mat) <- colnames(adj_mat) <- attr(adj_mat, "region.id") 
 
# Fit SpaMixed 
fit <- fit_SpaMixed(
   phy_obj = physeq_example,
   exposure = "PM2.5",
   covariates = c("age", "gender"),
   region_var = "zip_code",
   region_adj = adj_mat,
   micro_adj = adj_micro,
   scale_vars = c("PM2.5", "age"),
   factor_vars = c("gender"),
   lfdr_thres = 0.2
 )


```
