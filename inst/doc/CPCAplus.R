params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
library(multivarious)
library(ggplot2)

## ----motivation_example-------------------------------------------------------
set.seed(123)
n_samples <- 100
n_features <- 50

# Create background data (e.g., healthy controls)
# Main variation is in features 1-10
background <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
background[, 1:10] <- background[, 1:10] * 3  # Strong common variation

# Create foreground data (e.g., patients)
# Has the same common variation PLUS disease-specific signal in features 20-25
foreground <- background[1:60, ]  # Start with same structure
foreground[, 20:25] <- foreground[, 20:25] + matrix(rnorm(60 * 6, sd = 2), 60, 6)

# Standard PCA on combined data
all_data <- rbind(background, foreground)
regular_pca <- pca(all_data, ncomp = 2)

# Contrastive PCA
cpca_result <- cPCAplus(X_f = foreground, X_b = background, ncomp = 2)

# Compare what each method finds
loadings_df <- rbind(
  data.frame(
    feature = factor(1:30),
    value   = abs(regular_pca$v[1:30, 1]),
    method  = "Standard PCA"
  ),
  data.frame(
    feature = factor(1:30),
    value   = abs(cpca_result$v[1:30, 1]),
    method  = "Contrastive PCA"
  )
)

ggplot(loadings_df, aes(x = feature, y = value)) +
  geom_col(fill = "#1f78b4") +
  facet_wrap(~method, nrow = 1) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Feature",
    y = "|Loading|",
    title = "Top loading coefficients for PC1"
  )

## ----basic_usage--------------------------------------------------------------
# Basic usage
cpca_fit <- cPCAplus(
  X_f = foreground,  # Your group of interest (foreground)
  X_b = background,  # Your reference group (background)
  ncomp = 5          # Number of components to extract
)

# The result is a bi_projector object with familiar methods
print(cpca_fit)

# Project new data
new_samples <- matrix(rnorm(10 * n_features), 10, n_features)
new_scores <- project(cpca_fit, new_samples)

# Reconstruct using top components
reconstructed <- reconstruct(cpca_fit, comp = 1:2)

## ----understanding_output-----------------------------------------------------
# Which features contribute most to the first contrastive component?
top_features <- order(abs(cpca_fit$v[, 1]), decreasing = TRUE)[1:10]
print(paste("Top contributing features:", paste(top_features, collapse = ", ")))

# How much more variable is each component in foreground vs background?
print(paste("Variance ratios:", paste(round(cpca_fit$values[1:3], 2), collapse = ", ")))

## ----biomedical_example, eval=FALSE-------------------------------------------
# # Identify disease-specific patterns
# tumor_cpca <- cPCAplus(
#   X_f = tumor_samples,
#   X_b = healthy_tissue,
#   ncomp = 10
# )

## ----technical_example, eval=FALSE--------------------------------------------
# # Use technical replicates as background to find biological signal
# bio_cpca <- cPCAplus(
#   X_f = biological_samples,
#   X_b = technical_replicates,
#   ncomp = 5
# )

## ----time_example, eval=FALSE-------------------------------------------------
# # Find patterns specific to treatment timepoint
# treatment_cpca <- cPCAplus(
#   X_f = after_treatment,
#   X_b = before_treatment,
#   ncomp = 5
# )

## ----high_dim-----------------------------------------------------------------
# Create high-dimensional example
n_f <- 50; n_b <- 80; p <- 1000
X_background_hd <- matrix(rnorm(n_b * p), n_b, p)
X_foreground_hd <- X_background_hd[1:n_f, ] + 
                   matrix(c(rnorm(n_f * 20, sd = 2), rep(0, n_f * (p-20))), n_f, p)

# Use sample-space strategy for efficiency
cpca_hd <- cPCAplus(X_f = X_foreground_hd, X_b = X_background_hd,
                    ncomp = 5, strategy = "sample")

## ----regularization-----------------------------------------------------------
# Small background sample size can lead to instability
small_background <- matrix(rnorm(20 * 100), 20, 100)
small_foreground <- matrix(rnorm(30 * 100), 30, 100)

# Add regularization
cpca_regularized <- cPCAplus(
  X_f = small_foreground,
  X_b = small_background,
  ncomp = 5,
  lambda = 0.1  # Regularization parameter for background covariance
)

