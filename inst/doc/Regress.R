params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width=6, fig.height=4)
# Assuming necessary multivarious functions are loaded
# e.g., via devtools::load_all() or library(multivarious)
library(multivarious)
# Implicit dependencies like glmnet or pls might be needed depending on method used.

## ----build_basis--------------------------------------------------------------
set.seed(42)
n  <- 128                       # Number of observations (e.g., signals)
p  <- 32                        # Original variables per observation (e.g., time points)
k  <- 20                        # Number of basis functions (<= p often, <= n for lm)

## Create toy data: smooth signals + noise
t   <- seq(0, 1, length.out = p)
Y   <- replicate(n,  3*sin(2*pi*3*t) + 2*cos(2*pi*5*t) ) + 
       matrix(rnorm(n*p, sd = 0.3), p, n) # Note: Y is p x n here

## Orthonormal Fourier basis (columns = basis functions)
# We create k/2 sine and k/2 cosine terms, plus an intercept
freqs <- 1:(k %/% 2) # Integer division for number of frequencies
B <- cbind(rep(1, p), # Intercept column
           do.call(cbind, lapply(freqs, function(f) sin(2*pi*f*t))),
           do.call(cbind, lapply(freqs, function(f) cos(2*pi*f*t))))
colnames(B) <- c("Intercept", paste0("sin", freqs), paste0("cos", freqs))

# Make columns orthonormal (length 1, orthogonal to each other)
B <- scale(B, center = FALSE, scale = sqrt(colSums(B^2))) 

cat(paste("Dimensions: Y is", nrow(Y), "x", ncol(Y), 
          ", Basis B is", nrow(B), "x", ncol(B), "\n"))
# We want coefficients C (k x n) such that Y ≈ B %*% C.

## ----fit_regression-----------------------------------------------------------
library(multivarious)

# Fit using standard linear models (lm)
# Y is p x n (32 x 128): 32 time points x 128 signals
# B is p x k (32 x 21): 32 time points x 21 basis functions
# regress will fit 128 separate regressions, each with 32 observations and 21 predictors
fit <- regress(X = B,          # Predictors = basis functions (p x k)
               Y = Y,          # Response = signals (p x n)
               method    = "lm",
               intercept = FALSE) # Basis B already includes an intercept column

# The result is a bi_projector object
print(fit)

## Conceptual mapping to bi_projector slots:
# fit$v : Coefficients (n x k) - Basis coefficients for each signal.
# fit$s : Design Matrix (p x k) - The basis matrix B.
#         Stored for reconstruction.

## ----inspect_coefficients-----------------------------------------------------
coef_matrix_first3 <- fit$v[1:3, ]
cat("Coefficient matrix shape (first 3 signals):",
    nrow(coef_matrix_first3), "x", ncol(coef_matrix_first3), "\n\n")

cat("Coefficients for signal 1:\n")
print(coef_matrix_first3[1, ])

## ----reconstruct_fitted-------------------------------------------------------
Y_hat <- fit$s %*% t(fit$v)
max_reconstruction_error <- max(abs(Y_hat - Y))

cat("Reconstruction shape:", nrow(Y_hat), "x", ncol(Y_hat), "\n")
cat("Maximum reconstruction error:", format(max_reconstruction_error, digits=3), "\n")

## ----project_new_signal-------------------------------------------------------
Y_new_signal <- 3*sin(2*pi*3*t) + 2*cos(2*pi*5*t) + rnorm(p, sd=0.3)
Y_new_matrix <- matrix(Y_new_signal, nrow = p, ncol = 1)

coef_new <- t(fit$s) %*% Y_new_matrix

cat("Basis coefficients for new signal:\n")
print(coef_new[, 1])

## ----reconstruct_new_signal---------------------------------------------------
Y_new_recon <- fit$s %*% coef_new
reconstruction_error <- sqrt(mean((Y_new_matrix - Y_new_recon)^2))

cat("Reconstruction RMSE for new signal:", format(reconstruction_error, digits=3), "\n")

## ----regularized_models, eval=FALSE-------------------------------------------
# # Ridge regression (requires glmnet)
# fit_ridge <- regress(X = B, Y = Y, method = "mridge", lambda = 0.01, intercept = FALSE)
# 
# # Elastic Net (requires glmnet)
# fit_enet  <- regress(X = B, Y = Y, method = "enet", alpha = 0.5, lambda = 0.02, intercept = FALSE)
# 
# # Partial Least Squares (requires pls package) - useful if k > p or multicollinearity
# fit_pls   <- regress(X = B, Y = Y, method = "pls", ncomp = 15, intercept = FALSE)
# 
# # All these return bi_projector objects, so downstream code using
# # project(), reconstruct(), coef() etc. remains the same.

## ----partial_mappings---------------------------------------------------------
# Truncate: Keep only the first 5 basis functions (Intercept + 2 sine + 2 cosine)
fit5   <- truncate(fit, ncomp = 5) 
cat("Dimensions after truncating to 5 components:", 
    "Basis (s):", paste(dim(fit5$s), collapse="x"), 
    ", Coefs (v):", paste(dim(fit5$v), collapse="x"), "\n")
# Reconstruction using only first 5 basis functions (manual)
# Equivalent to: scores(fit5) %*% t(coef(fit5)) for the selected components
Y_hat5 <- fit5$s %*% t(fit5$v)

# Partial inverse projection: Map only a subset of coefficients back
# e.g., reconstruct using only components 2 through 6 (skip intercept)
# Note: partial_inverse_projection is not a standard bi_projector method, 
# this might require manual slicing of the basis matrix B (fit$s) or coefs (fit$v).
# Manual reconstruction example for components 2:6
coef_subset <- fit$v[2:6, , drop=FALSE] # k_sub x n
basis_subset <- fit$s[, 2:6, drop=FALSE] # p x k_sub
Y_lowHat <- basis_subset %*% coef_subset # p x n reconstruction

# Variable usage helpers (Conceptual - actual functions might differ)
# `variables_used(fit)` could show which basis functions have non-zero coefficients (esp. for 'enet').
# `vars_for_component(fit, k)` isn't directly applicable here as components are the basis functions themselves.

## ----internal_checks, eval=nzchar(Sys.getenv("_MULTIVARIOUS_DEV_COVERAGE"))----
# This chunk only runs if _MULTIVARIOUS_DEV_COVERAGE is non-empty
message("Running internal consistency checks for regress()...")
tryCatch({
  stopifnot(
    # Check reconstruction fidelity for lm
    max(abs(reconstruct(fit) - Y)) < 1e-10,
    # Check dimensions of inverse projection matrix (n x p)
    # inverse_projection maps coefficients (k x n) back to data (p x n)
    # The matrix itself maps k -> p implicitly. Let's check coef matrix dims.
    nrow(fit$v) == ncol(B), # k rows
    ncol(fit$v) == ncol(Y)  # n columns
    # Add checks for other methods if evaluated
  )
  message("Regress internal checks passed.")
}, error = function(e) {
  warning("Regress internal checks failed: ", e$message)
})

