params <-
list(family = "red")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 4
)
library(multivarious)
library(dplyr)
library(ggplot2)

## ----data_partial_proj--------------------------------------------------------
set.seed(1)
n  <- 100
p  <- 8
X  <- matrix(rnorm(n * p), n, p)

# Fit a centred 3-component PCA (via SVD)
# Manually center the data and create fitted preprocessor
Xc <- scale(X, center = TRUE, scale = FALSE)
svd_res <- svd(Xc, nu = 0, nv = 3)

# Create a fitted centering preprocessor
preproc_fitted <- fit(center(), X)

pca <- bi_projector(
  v     = svd_res$v,
  s     = Xc %*% svd_res$v,
  sdev  = svd_res$d[1:3] / sqrt(n-1), # Correct scaling for sdev
  preproc = preproc_fitted
)

## ----project_full-------------------------------------------------------------
scores_full <- project(pca, X)        # n × 3
head(round(scores_full, 2))

## ----project_partial----------------------------------------------------------
X_miss      <- X[, 1:6]               # keep only first 6 columns
col_subset  <- 1:6                    # their positions in the **original** X

scores_part <- partial_project(pca, X_miss, colind = col_subset)

# How close are the results?
plot_df <- tibble(
  full = scores_full[,1],
  part = scores_part[,1]
)

ggplot(plot_df, aes(full, part)) +
  geom_point() +
  geom_abline(col = "red") +
  coord_equal() +
  labs(title = "Component 1: full vs. partial projection") +
  theme_minimal()

## ----partial_projector_cache--------------------------------------------------
# Assuming partial_projector is available
pca_1to6 <- partial_projector(pca, 1:6)   # keeps a reference + cache

# project 1000 new observations that only have the first 6 vars
new_batch <- matrix(rnorm(1000 * 6), 1000, 6)
scores_fast <- project(pca_1to6, new_batch)
dim(scores_fast)   # 1000 × 3

## ----multiblock_example-------------------------------------------------------
# Create a multiblock projector from our PCA
# Suppose columns 1-4 are "Block A" (block 1) and columns 5-8 are "Block B" (block 2)
block_indices <- list(1:4, 5:8)

mb <- multiblock_projector(
  v = pca$v,
  preproc = pca$preproc,
  block_indices = block_indices
)

# Now we can project using only Block 2's data (columns 5-8)
X_block2 <- X[, 5:8]
scores_block2 <- project_block(mb, X_block2, block = 2)

# Compare to full projection
head(round(cbind(full = scores_full[,1], block2 = scores_block2[,1]), 2))

## ----roi_project--------------------------------------------------------------
roi_cols   <- 1:5                 # pretend these are the ROI voxels
X_roi      <- X[, roi_cols]       # same matrix from Section 2

roi_scores <- partial_project(pca, X_roi, colind = roi_cols)

# Compare component 1 from full vs ROI
df_roi <- tibble(
  full = scores_full[,1],
  roi  = roi_scores[,1]
)

ggplot(df_roi, aes(full, roi)) +
  geom_point(alpha = .6) +
  geom_abline(col = "red") +
  coord_equal() +
  labs(title = "Component 1 scores: full data vs ROI") +
  theme_minimal()


## ----block_single_subject-----------------------------------------------------
# Get scores for observation 1 using only Block 1 variables (columns 1-4)
subject1_block1 <- project_block(mb, X[1, 1:4, drop = FALSE], block = 1)

# Get scores for the same observation using only Block 2 variables (columns 5-8)
subject1_block2 <- project_block(mb, X[1, 5:8, drop = FALSE], block = 2)

# Compare: do both blocks tell the same story about this observation?
cat("Subject 1 scores from Block 1:", round(subject1_block1, 2), "\n")
cat("Subject 1 scores from Block 2:", round(subject1_block2, 2), "\n")
cat("Subject 1 scores from full data:", round(scores_full[1,], 2), "\n")

## ----session-info-extra-------------------------------------------------------
sessionInfo()

