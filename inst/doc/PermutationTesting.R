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

## ----mixed_effect_perm_example------------------------------------------------
set.seed(11)

design <- expand.grid(
  subject = factor(seq_len(8)),
  level = factor(c("low", "mid", "high"), levels = c("low", "mid", "high")),
  KEEP.OUT.ATTRS = FALSE
)
design$group <- factor(rep(c("A", "B"), each = 12))

level_num <- c(low = -1, mid = 0, high = 1)[as.character(design$level)]
group_num <- ifelse(design$group == "B", 1, 0)
subj_idx <- as.integer(design$subject)
b0 <- rnorm(8, sd = 0.5)

Y <- cbind(
  b0[subj_idx] + level_num + rnorm(nrow(design), sd = 0.15),
  group_num + rnorm(nrow(design), sd = 0.15),
  level_num * group_num + rnorm(nrow(design), sd = 0.15),
  rnorm(nrow(design), sd = 0.15)
)

fit_mixed <- mixed_regress(
  Y,
  design = design,
  fixed = ~ group * level,
  random = ~ 1 | subject,
  basis = shared_pca(3),
  preproc = pass()
)

E_int <- effect(fit_mixed, "group:level")
pt_int <- perm_test(E_int, nperm = 49, alpha = 0.10)

pt_int$component_results

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

