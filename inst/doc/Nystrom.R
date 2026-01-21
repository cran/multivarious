params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width=6, fig.height=4)
library(multivarious)
library(stats)

## ----quick_start--------------------------------------------------------------
set.seed(42)
N <- 1000
p <- 50
X <- matrix(rnorm(N * p), N, p)

# Nyström approximation with linear kernel (default)
# Uses only 100 landmarks instead of all 1000 points
ny_fit <- nystrom_approx(

  X,
  ncomp = 10,
  nlandmarks = 100,
  preproc = center()
)

print(ny_fit)

# Standard bi_projector interface
head(scores(ny_fit))

## ----custom_kernel------------------------------------------------------------
# RBF (Gaussian) kernel
rbf_kernel <- function(X, Y = NULL, sigma = 1) {
  if (is.null(Y)) Y <- X

  sumX2 <- rowSums(X^2)
  sumY2 <- rowSums(Y^2)
  sqdist <- outer(sumX2, sumY2, `+`) - 2 * tcrossprod(X, Y)
  sqdist[sqdist < 0] <- 0

  exp(-sqdist / (2 * sigma^2))
}

ny_rbf <- nystrom_approx(
  X,
  kernel_func = rbf_kernel,
  ncomp = 10,
  nlandmarks = 100
)

## ----projection---------------------------------------------------------------
X_new <- matrix(rnorm(50 * p), 50, p)
new_scores <- project(ny_fit, X_new)
dim(new_scores)

## ----double_nystrom-----------------------------------------------------------
# Standard method
system.time(

  ny_standard <- nystrom_approx(X, ncomp = 5, nlandmarks = 200, method = "standard")
)

# Double Nyström (faster with intermediate rank l)
system.time(
  ny_double <- nystrom_approx(X, ncomp = 5, nlandmarks = 200, method = "double", l = 50)
)

