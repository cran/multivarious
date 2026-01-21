params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(multivarious)
# library(future) # Load if using parallel = TRUE
# library(MASS) # Load if using default fit_fun for discriminant_projector

## ----basic_workflow, message=FALSE, warning=FALSE-----------------------------
data(iris)
X_iris <- as.matrix(iris[, 1:4])

mod_pca <- pca(X_iris, ncomp = 4, preproc = center())

## ----perm_test_call-----------------------------------------------------------
set.seed(1)
pt_pca <- perm_test(mod_pca,
                    X = X_iris,
                    nperm = 199,
                    comps = 3,
                    parallel = FALSE)

## ----inspect_results----------------------------------------------------------
print(pt_pca$component_results)

## ----custom_measure-----------------------------------------------------------
my_pca_stat <- function(model_perm, comp_idx, ...) {
  # Only compute the joint statistic when testing component 2

  if (comp_idx == 2 && length(model_perm$sdev) >= 2) {
    sum(model_perm$sdev[1:2]^2) / sum(model_perm$sdev^2)
  } else if (comp_idx == 1) {
    model_perm$sdev[1]^2 / sum(model_perm$sdev^2)
  } else {
    NA_real_
  }
}

# Illustrative call (using default measure here for simplicity)
pt_pca_custom <- perm_test(mod_pca, X = X_iris, nperm = 50, comps = 2,
                           parallel = FALSE)
print(pt_pca_custom$component_results)

## ----parallel_example, eval=FALSE---------------------------------------------
# library(future)
# plan(multisession, workers = 4)
# 
# pt_pca_parallel <- perm_test(mod_pca, X = X_iris,
#                              nperm = 999,
#                              comps = 3,
#                              parallel = TRUE)
# 
# plan(sequential)

## ----internal_checks, eval=nzchar(Sys.getenv("_MULTIVARIOUS_DEV_COVERAGE")), include=FALSE----
# CI sanity check: verify perm_test returns expected structure
set.seed(42)
mtcars_mat <- as.matrix(scale(mtcars))
pca_test_mod <- pca(mtcars_mat, ncomp = 3)
pt_check <- perm_test(pca_test_mod, mtcars_mat, nperm = 19, comps = 2, parallel = FALSE)
stopifnot(

  nrow(pt_check$component_results) == 2,
 !is.na(pt_check$component_results$pval[1])
)

