params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
library(multivarious)  # Ensure pca(), svd_wrapper(), bi_projector, etc. are available
library(ggplot2)       # for plots below

## ----svd_wrapper_example------------------------------------------------------
set.seed(1)
X <- matrix(rnorm(35*10), 35, 10)   # 35 obs × 10 vars

sv_fast <- svd_wrapper(X, ncomp = 5, preproc = center(), method = "fast")

# irlba backend (if installed) gives identical results
sv_irlba <- if (requireNamespace("irlba", quietly = TRUE)) {
  suppressWarnings(svd_wrapper(X, ncomp = 5, preproc = center(), method = "irlba"))
}

# Same downstream code works for both objects:
head(scores(sv_fast)) # 35 × 5

if (!is.null(sv_irlba)) {
  all.equal(scores(sv_fast), scores(sv_irlba))
}

## ----pca_example--------------------------------------------------------------
data(iris)
X_iris <- as.matrix(iris[, 1:4])

pca_fit <- pca(X_iris, ncomp = 4)    # defaults to method = "fast", preproc=center()
print(pca_fit)

## ----pca_screeplot------------------------------------------------------------
screeplot(pca_fit, type = "lines", main = "Iris PCA – scree plot")

## ----pca_biplot, warning=FALSE, message=FALSE---------------------------------
# Requires ggrepel for repulsion, but works without it
biplot(pca_fit, repel_points = TRUE, repel_vars = TRUE)

## ----biprojector_methods------------------------------------------------------
# rank-2 reconstruction of the iris data
Xhat2 <- reconstruct(pca_fit, comp = 1:2)
print(paste("MSE (rank 2):", round(mean((X_iris - Xhat2)^2), 4))) # MSE ~ 0.076

# drop to 2 PCs everywhere
pca2 <- truncate(pca_fit, 2)
shape(pca2)            # 4 vars × 2 comps

## ----code_coverage------------------------------------------------------------
# std scores
head(std_scores(svd_wrapper(X, ncomp = 3))) # Use the earlier X data

# tiny permutation test (10 perms; obviously too few for science)
# This requires perm_test.pca method
# Make sure X_iris is centered if perm_test needs centered data
perm_res <- perm_test(pca_fit, X_iris, nperm = 10, comps = 2)
print(perm_res$component_results)

# quick varimax rotation
if (requireNamespace("GPArotation", quietly = TRUE)) {
  pca_rotated <- rotate(pca_fit, ncomp = 3, type = "varimax")
  print(pca_rotated)
} else {
  cat("GPArotation not installed, skipping rotation example.\n")
}

## ----session_info_svd---------------------------------------------------------
sessionInfo()

