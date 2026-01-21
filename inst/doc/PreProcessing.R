params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width=6, fig.height=4)
library(multivarious)
library(dplyr) # Needed for %>% and tibble manipulation
library(tibble)
library(ggplot2)

## ----setup_data_preproc-------------------------------------------------------
set.seed(0)
X <- matrix(rnorm(10*4), 10, 4)

pp_pass <- fit(pass(), X)        # == do nothing
Xp_pass <- transform(pp_pass, X) # applies nothing, just copies X
all.equal(Xp_pass, X)            # TRUE

## ----standardize_example------------------------------------------------------
# Fit the preprocessor (calculates means & SDs from X) and transform
pp_std <- fit(standardize(), X)
Xs     <- transform(pp_std, X)

# Check results
all(abs(colMeans(Xs)) < 1e-12)   # TRUE: data is centered
round(apply(Xs, 2, sd), 6)       # ~1: data is scaled

# Check back-transform
all.equal(inverse_transform(pp_std, Xs), X) # TRUE

## ----partial_transform--------------------------------------------------------
X_cols24 <- X[, c(2,4), drop=FALSE] # Keep as matrix

# Apply the *already fitted* standardizer using only columns 2 & 4
Xs_cols24 <- transform(pp_std, X_cols24, colind = c(2,4))

# Compare original columns 2, 4 with their transformed versions
head(cbind(X_cols24, Xs_cols24))

# Back-transform works too
X_rev_cols24 <- inverse_transform(pp_std, Xs_cols24, colind = c(2,4))
all.equal(X_rev_cols24, X_cols24) # TRUE

## ----pipe_example-------------------------------------------------------------
# Define a pipeline: center, then scale to unit variance
# Fit the pipeline to the data
pp_pipe <- fit(standardize(), X)

# Apply the pipeline
Xp_pipe <- transform(pp_pipe, X)

## ----plot_pipeline------------------------------------------------------------
# Compare first column before and after pipeline
df_pipe <- tibble(raw = X[,1],   processed = Xp_pipe[,1])

ggplot(df_pipe) +
  geom_density(aes(raw), colour = "red", linewidth = 1) +
  geom_density(aes(processed), colour = "blue", linewidth = 1) +
  ggtitle("Column 1 Density: Before (red) and After (blue) Pipeline") +
  theme_minimal()

## ----concat_example-----------------------------------------------------------
# Two fake blocks with distinct scales
X1 <- matrix(rnorm(10*5 , 10 , 5), 10, 5)   # block 1: high mean
X2 <- matrix(rnorm(10*7 ,  2 , 7), 10, 7)   # block 2: low mean

# Fit separate preprocessors for each block
p1 <- fit(center(), X1)
p2 <- fit(standardize(), X2)

# Transform each block
X1p <- transform(p1, X1)
X2p <- transform(p2, X2)

# Concatenate the *fitted* preprocessors
block_indices_list = list(1:5, 6:12)
pp_concat <- concat_pre_processors(
  list(p1, p2),
  block_indices = block_indices_list
)

# Apply the concatenated preprocessor to the combined data
X_combined <- cbind(X1, X2)
X_combined_p <- transform(pp_concat, X_combined)

# Check means (block 1 only centered, block 2 standardized)
round(colMeans(X_combined_p), 2)

# Need only block 1 processed later? Use colind with global indices
X1_later_p <- transform(pp_concat, X1, colind = block_indices_list[[1]])
all.equal(X1_later_p, X1p) # TRUE

# Need block 2 processed?
X2_later_p <- transform(pp_concat, X2, colind = block_indices_list[[2]])
all.equal(X2_later_p, X2p) # TRUE

## ----concat_reversibility-----------------------------------------------------
back_combined <- inverse_transform(pp_concat, X_combined_p)

# Compare first few rows/cols of original vs round-trip
knitr::kable(
  head(cbind(orig = X_combined[, 1:6], recon = back_combined[, 1:6]), 3),
  digits = 2,
  caption = "First 3 rows, columns 1-6: Original vs Reconstructed"
)

all.equal(X_combined, back_combined) # TRUE

## ----session_info_preproc-----------------------------------------------------
sessionInfo()

