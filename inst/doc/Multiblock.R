params <-
list(family = "red")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 7,
  fig.height = 4
)
library(dplyr)
library(multivarious)
# Assuming necessary multiblock functions are loaded, e.g., via devtools::load_all()

## ----data_multiblock----------------------------------------------------------
set.seed(1)
n  <- 100
pA <- 7; pB <- 5                    # two blocks, different widths

XA <- matrix(rnorm(n * pA), n, pA)
XB <- matrix(rnorm(n * pB), n, pB)
X  <- cbind(XA, XB)                 # global data matrix
blk_idx <- list(A = 1:pA, B = (pA + 1):(pA + pB)) # Named list is good practice

## ----build_multiblock---------------------------------------------------------
# 2-component centred PCA (using base SVD for brevity)
preproc_fitted <- fit(center(), X)
Xc        <- transform(preproc_fitted, X)          # Centered data
svd_res   <- svd(Xc, nu = 0, nv = 2)               # only V (loadings)
mb        <- multiblock_projector(
  v             = svd_res$v,                       # p × k loadings
  preproc       = preproc_fitted,                  # remembers centering
  block_indices = blk_idx
)

print(mb)

## ----project_multiblock_all---------------------------------------------------
scores_all <- project(mb, X)                       # n × 2
head(round(scores_all, 3))

## ----project_multiblock_block-------------------------------------------------
# Project using only data from block A (requires original columns)
scores_A <- project_block(mb, XA, block = 1)       
# Project using only data from block B
scores_B <- project_block(mb, XB, block = 2)       

cor(scores_all[,1], scores_A[,1])                  # high (they coincide)

## ----project_multiblock_partial-----------------------------------------------
# Get the global indices for the first 3 columns of block B
sel_cols_global <- blk_idx[["B"]][1:3]
# Extract the corresponding data columns from the full matrix or block B
part_XB_data  <- X[, sel_cols_global, drop = FALSE] # Data must match global indices

scores_part <- partial_project(mb, part_XB_data,
                               colind = sel_cols_global)  # Use global indices
head(round(scores_part, 3))

## ----build_biprojector--------------------------------------------------------
bi <- multiblock_biprojector(
  v             = svd_res$v,
  s             = Xc %*% svd_res$v,    # Calculate scores: Xc %*% V
  sdev          = svd_res$d[1:2] / sqrt(n-1), # SVD d are related to sdev
  preproc       = preproc_fitted,
  block_indices = blk_idx
)
print(bi)

## ----perm_test_multiblock-----------------------------------------------------
# Quick permutation test (use more permutations for real analyses)
# use_rspectra=FALSE needed for this 2-block example; larger problems can use TRUE
perm_res <- perm_test(bi, Xlist = list(A = XA, B = XB), nperm = 99, use_rspectra = FALSE)
print(perm_res$component_results)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

