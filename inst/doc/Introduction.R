## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

# Load the multivarious package
library(multivarious)

# Create a synthetic dataset
set.seed(42)
X <- matrix(rnorm(200), 10, 20)

# Perform SVD on the dataset
svdfit <- svd(X)

# Create a bi_projector object
p <- bi_projector(svdfit$v, s = svdfit$u %*% diag(svdfit$d), sdev = svdfit$d)

# Generate new data to project onto the same subspace as the original data
new_data <- matrix(rnorm(5 * 20), 5, 20)

projected_data <- project(p, new_data)
print(projected_data)


## -----------------------------------------------------------------------------
# Load iris dataset and select the first four columns
data(iris)
X <- iris[, 1:4]

# Compute SVD using the base method and 3 components
fit <- svd_wrapper(X, ncomp = 3, preproc = center(), method = "base")


## -----------------------------------------------------------------------------
# Define new_data
new_data <- rnorm(nrow(iris))

# Project the new variables into the subspace
projected_vars <- project_vars(fit, new_data)


