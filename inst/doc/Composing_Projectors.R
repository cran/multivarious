params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width=6, fig.height=4)
# Assuming necessary multivarious functions are loaded
# e.g., via devtools::load_all() or library(multivarious)
library(multivarious)
library(tibble) # For summary output

## ----quick_start--------------------------------------------------------------
set.seed(1)

X  <- matrix(rnorm(30*15), 30, 15)   # raw data, 30 samples, 15 variables
p1 <- pca(X, ncomp = 8)              # first reduction: 15 -> 8 components
p2 <- pca(scores(p1), ncomp = 7)     # second reduction: 8 -> 4 components

# Compose the two projectors
pipe <- compose_partial_projector(
           first  = p1,
           second = p2)

print(pipe)

# Project original data through the entire pipeline
S <- project(pipe, X)          # 30 × 4 scores – as if the two steps were one
dim(S)

# Get a summary of the pipeline stages
summary(pipe)

## ----partial_projection_examples----------------------------------------------
# Example 1: Use only variables 1:5 for the *first* PCA stage.
# The second PCA stage receives the full 8 components from the (partial) first stage.
S15 <- partial_project(pipe, X[, 1:5, drop=FALSE], colind = 1:5)
cat("Dimensions after partial projection (cols 1:5 in first stage):", dim(S15), "\n")

# Example 2: Multi-stage pipeline (conceptual)
# Imagine a 3-stage pipeline: wavelets -> PCA (block1) -> PCA (global)
# pipe2 <- wavelet_projector(...) %>>% 
#          pca(..., ncomp = 10)   %>>% 
#          pca(..., ncomp = 3)

# To focus on coefficients 12:20 *after* the wavelet step (i.e., input to stage 2):
# S_sel <- partial_project(pipe2, X, # Assuming X is appropriate input for wavelets
#                          colind = list(NULL, 12:20, NULL))
# Note: The indices in the list always refer to the dimensions *entering* that specific stage.

## ----reconstruction-----------------------------------------------------------
# Reconstruct original data from the final scores 'S'
X_hat <- reconstruct(pipe, S)
cat("Dimensions of reconstructed data:", dim(X_hat), "\n")

# Check reconstruction accuracy
# Note: Since the pipeline involves dimensionality reduction (15 -> 8 -> 4),
# reconstruction will not be exact. The error reflects the information lost.
max_reconstruction_error <- max(abs(X - X_hat))
cat("Maximum absolute reconstruction error:", format(max_reconstruction_error, digits=3), "\n")
# stopifnot(max_reconstruction_error < 1e-5) # Removed: This check is too strict for lossy reconstruction

# Get the overall coefficient matrix (p_orig x q_final)
V <- coef(pipe)
cat("Dimensions of overall coefficient matrix:", dim(V), "\n")

# Get the overall pseudo-inverse matrix (q_final x p_orig)
Vplus <- inverse_projection(pipe)
cat("Dimensions of overall inverse projection matrix:", dim(Vplus), "\n")

## ----helper_pipe, eval=FALSE--------------------------------------------------
# # pipe3 <- pca1 %>>% pca2 %>>% pca3

