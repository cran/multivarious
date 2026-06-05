params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(multivarious)

## ----pca_example--------------------------------------------------------------
# Load iris dataset and select numeric columns
data(iris)
X <- as.matrix(iris[, 1:4])

# 1. Define a pre-processor (center the data)
preproc <- center()

# 2. Fit PCA using svd_wrapper, keeping 3 components
#    The pre-processor is applied internally.
fit <- pca(X, ncomp = 3, preproc = preproc)

# The result 'fit' is a bi_projector
print(fit)

# 3. Access results
iris_scores <- scores(fit) # Scores of the centered training data (150 x 3)
iris_loadings <- loadings(fit) # Loadings (4 x 3)
cat("\nDimensions of Scores:", dim(iris_scores), "\n")
cat("Dimensions of Loadings:", dim(iris_loadings), "\n")

# 4. Project new data
# Create some new iris-like samples (5 samples, 4 variables)
set.seed(123)
new_iris_data <- matrix(rnorm(5 * 4, mean = colMeans(X), sd = apply(X, 2, sd)), 
                        nrow = 5, byrow = TRUE)

# Project the new data into the PCA space defined by 'fit'
# Pre-processing (centering using training data means) is applied automatically.
projected_new_scores <- project(fit, new_iris_data)
cat("\nDimensions of Projected New Data Scores:", dim(projected_new_scores), "\n")
print(head(projected_new_scores))

# 5. Reconstruct approximated original data from scores
# Reconstruct the first few original samples
reconstructed_X_approx <- reconstruct(fit, comp=1:3) # uses scores(fit) by default
cat("\nReconstructed Approximation of Original Data (first 5 rows):\n")
print(head(reconstructed_X_approx))

print(head(X)) # Original data for comparison


## ----mixed_effect_example-----------------------------------------------------
set.seed(99)

design_m <- expand.grid(
  subject = factor(seq_len(6)),
  level = factor(c("low", "mid", "high"), levels = c("low", "mid", "high")),
  KEEP.OUT.ATTRS = FALSE
)
design_m$group <- factor(rep(c("A", "B"), each = 9))

level_num <- c(low = -1, mid = 0, high = 1)[as.character(design_m$level)]
group_num <- ifelse(design_m$group == "B", 1, 0)
subj_idx <- as.integer(design_m$subject)
b0 <- rnorm(6, sd = 0.5)

Y_m <- cbind(
  b0[subj_idx] + level_num + rnorm(nrow(design_m), sd = 0.15),
  group_num + rnorm(nrow(design_m), sd = 0.15),
  level_num * group_num + rnorm(nrow(design_m), sd = 0.15),
  rnorm(nrow(design_m), sd = 0.15)
)

fit_m <- mixed_regress(
  Y_m,
  design = design_m,
  fixed = ~ group * level,
  random = ~ 1 | subject,
  basis = shared_pca(3),
  preproc = pass()
)

E_gl <- effect(fit_m, "group:level")
pt_gl <- perm_test(E_gl, nperm = 19, alpha = 0.10)

print(E_gl)
pt_gl$component_results

## ----project_vars_example-----------------------------------------------------
# Using the 'fit' object from the PCA example above

# Create a new variable (column) with the same number of samples as original data
set.seed(456)
new_variable <- rnorm(nrow(X))

# Project this new variable into the component space defined by the PCA scores (fit$s)
# Result shows how the new variable relates to the principal components.
projected_variable_loadings <- project_vars(fit, new_variable)
cat("\nProjection of new variable onto components:", projected_variable_loadings, "\n")

