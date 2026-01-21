params <-
list(family = "red")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 5
)
# Load necessary packages for the examples
# Ensure 'multiblock' package itself is loaded for testing vignettes,
# e.g., via devtools::load_all() or library(multiblock)
library(tibble)
library(dplyr)
library(stats) # For cancor
library(glmnet) # For glmnet example
library(multivarious)
# Helper k-fold function (replace with package internal if available)
kfold_split <- function(n, k = 5) {
  idx <- sample(rep(1:k, length.out = n))
  lapply(1:k, function(j) list(train = which(idx != j),
                               test  = which(idx == j)))
}

## ----data_cca-----------------------------------------------------------------
data(iris)
X <- as.matrix(iris[, 1:2])
Y <- as.matrix(iris[, 3:4])

# Show first few rows of combined data
head(cbind(X, Y))

## ----fit_cca_wrapper----------------------------------------------------------
fit_cca <- function(Xtrain, Ytrain, ncomp = 2, ...) {
  # Step 1: Define and fit preprocessors with training data
  # Use center()+z-score for both blocks (swap to center() only if preferred)
  preproc_x <- fit(colscale(center(), type = "z"), Xtrain)
  preproc_y <- fit(colscale(center(), type = "z"), Ytrain)

  # Step 2: Transform training data with the fitted preprocessors
  Xp <- transform(preproc_x, Xtrain)
  Yp <- transform(preproc_y, Ytrain)

  # Step 3: Run CCA on the preprocessed data (no additional centering here)
  cc <- stats::cancor(Xp, Yp, xcenter = FALSE, ycenter = FALSE)

  # Step 4: Store loadings and preprocessors in a cross_projector
  cp <- cross_projector(
    vx         = cc$xcoef[, 1:ncomp, drop = FALSE],
    vy         = cc$ycoef[, 1:ncomp, drop = FALSE],
    preproc_x  = preproc_x,
    preproc_y  = preproc_y,
    classes    = "cca_cross_projector"
  )
  attr(cp, "can_cor") <- cc$cor[1:ncomp]
  cp
}

# Fit the model
cp_cca <- fit_cca(X, Y)
print(cp_cca)
attr(cp_cca, "can_cor") # Show canonical correlations

## ----transfer_cca-------------------------------------------------------------
# Predict Y from X
Y_hat <- transfer(cp_cca, X, from = "X", to = "Y")
head(round(Y_hat, 2))

# Evaluate reconstruction quality (comparing original Y to predicted Y)
# Note: CCA maximizes correlation in the latent space, not prediction accuracy.
# With only 2 components and weak relationships between sepal and petal features,
# the R² may be modest. For prediction tasks, consider PLS instead of CCA.
measure_reconstruction_error(Y, Y_hat, metrics = c("rmse", "r2"))

## ----partial_project_cca------------------------------------------------------
new_x_partial   <- X[, 1, drop = FALSE] # One column matrix
col_index       <- 1                     # Its index in the original X block

# Project this partial data into the latent space
# Note: partial_project needs 'source' for cross_projector
scores_partial  <- partial_project(cp_cca, new_x_partial,
                                   colind = col_index,
                                   source = "X") # Specify source block
head(round(scores_partial, 3))

# We could then map these scores back to the Y-block space if needed,
# though interpretation might require care. Example using inverse_projection:
# Assuming inverse_projection method exists for cross_projector
# y_from_partial  <- scores_partial %*% inverse_projection(cp_cca, domain = "Y")
# head(round(y_from_partial, 2))

## ----cv_cca-------------------------------------------------------------------
set.seed(1)
folds <- kfold_split(nrow(X), k = 5) # Create 5 folds

# Note: Using cv_generic requires careful handling of two-block data.
# The following demonstrates the concept, assuming a cv method exists
# that can handle X, Y inputs directly or that cv_generic is adapted.
# For a real application, one might need helper functions within .fit_fun
# and .measure_fun to split/combine the data blocks within each fold.

# Placeholder call assuming an appropriate `cv` method or adapted `cv_generic`:
# cv_res <- cv( # Or cv_generic(...)
#   X, Y, # Hypothetical interface accepting two blocks
#   folds          = folds,
#   max_comp       = 2,
#   fit_fun        = fit_cca,
#   measure_fun    = measure_interblock_transfer_error,
#   measure_args   = list(metrics = c("x2y.rmse", "y2x.rmse"))
# )
# print(summary(cv_res)) # Would show average RMSE per component

# Example using cv_generic with wrappers (conceptual)
# Define wrappers carefully based on how data/folds are structured
# .fit_wrapper <- function(train_data_list, ncomp, ...) {
#    fit_cca(train_data_list$X, train_data_list$Y, ncomp=ncomp, ...)
# }
# .measure_wrapper <- function(model, test_data_list, metrics) {
#    measure_interblock_transfer_error(test_data_list$X, test_data_list$Y, model, metrics=metrics)
# }
# cv_generic_res <- cv_generic(
#    data = list(X=X, Y=Y), # Requires folds based on row indices of X & Y
#    folds = folds,
#    .fit_fun = .fit_wrapper,
#    fit_args = list(ncomp=2), # Pass ncomp here if fit_cca needs it
#    .measure_fun = .measure_wrapper,
#    measure_args = list(metrics=c("x2y.rmse", "y2x.rmse"))
# )
# print(summary(cv_generic_res)) # This would work if wrappers are correct

cat("Skipping CV execution in vignette for brevity/simplicity.\n")
cat("Users should adapt cv/cv_generic call based on package structure.\n")

## ----perm_test_cca------------------------------------------------------------
# Assuming perm_test.cross_projector method exists
# perm_res <- perm_test(
#   cp_cca, X, Y = Y, # Pass the cross_projector and original data blocks
#   nperm = 100, # Use more perms (e.g., 1000) in practice
#   alternative = "less", # Test if observed error is less than permuted error
#   # Need to specify the statistic measured by perm_test for cross_projector
#   # e.g., measure_fun = function(model, X, Y) {
#   #    Yh <- transfer(model, X, from="X", to="Y"); measure_reconstruction_error(Y, Yh, "rmse")$rmse
#   # }
# )
# print(perm_res)

cat("Skipping permutation test execution in vignette.\n")
cat("Requires perm_test method for cross_projector, specifying the test statistic.\n")

## ----setup_glmnet-------------------------------------------------------------
# Generate sample data
set.seed(123)
n_obs <- 100
n_pred <- 50
X_glm <- matrix(rnorm(n_obs * n_pred), n_obs, n_pred)
# True coefficients (sparse)
true_beta <- matrix(0, n_pred, 1)
true_beta[1:10, 1] <- runif(10, -1, 1)
# Response variable with noise
y_glm <- X_glm %*% true_beta + rnorm(n_obs, sd = 0.5)

# Fit glmnet (LASSO, alpha=1)
# Typically use cv.glmnet to find lambda, but using a fixed one for simplicity
glm_fit <- glmnet::glmnet(X_glm, y_glm, alpha = 1) # alpha=1 is LASSO

# Choose a lambda (e.g., one near the end of the path)
chosen_lambda <- glm_fit$lambda[length(glm_fit$lambda) * 0.8]

# Get coefficients for this lambda
beta_hat <- coef(glm_fit, s = chosen_lambda)
print(paste("Number of non-zero coefficients:", sum(beta_hat != 0)))

# Extract coefficients, excluding the intercept
v_glm <- beta_hat[-1, 1, drop = FALSE] # Drop intercept, ensure it's a matrix
dim(v_glm) # Should be n_pred x 1

## ----wrap_glmnet--------------------------------------------------------------
# Define preprocessor
# Define and fit preprocessor with training data
preproc_glm <- fit(colscale(center(), type = "z"), X_glm)

# Create the projector
proj_glm <- projector(
  v = v_glm,
  preproc = preproc_glm,
  classes = "glmnet_projector"
)

print(proj_glm)

## ----project_glmnet-----------------------------------------------------------
# Generate some new test data
X_glm_test <- matrix(rnorm(20 * n_pred), 20, n_pred)

# Project the test data
# project() handles applying the centering/scaling from preproc_glm
lasso_scores <- project(proj_glm, X_glm_test)

head(round(lasso_scores, 3))
dim(lasso_scores) # Should be n_test x 1 (since v_glm has 1 column)

# Compare with direct calculation using transform
# Apply the same preprocessing used within project()
X_glm_test_processed <- transform(preproc_glm, X_glm_test)
# Calculate scores directly using the processed data and coefficients
direct_scores <- X_glm_test_processed %*% v_glm
head(round(direct_scores, 3))

# Check they are close
all.equal(c(lasso_scores), c(direct_scores))

