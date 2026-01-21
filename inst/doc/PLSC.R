params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5)
library(multivarious)
library(ggplot2)
set.seed(1)

## ----simulate-----------------------------------------------------------------
n  <- 80   # subjects
pX <- 8    # brain/features
pY <- 5    # behavior
d  <- 2    # true latent dimensions

# orthonormal loadings
Vx_true <- qr.Q(qr(matrix(rnorm(pX * d), pX, d)))
Vy_true <- qr.Q(qr(matrix(rnorm(pY * d), pY, d)))

F_scores <- matrix(rnorm(n * d), n, d)         # latent scores
noise    <- 0.10

X <- F_scores %*% t(Vx_true) + noise * matrix(rnorm(n * pX), n, pX)
Y <- F_scores %*% t(Vy_true) + noise * matrix(rnorm(n * pY), n, pY)

## ----fit----------------------------------------------------------------------
fit_plsc <- plsc(X, Y, ncomp = 3,                   # request a few extra comps
                 preproc_x = standardize(),         # correlation-scale
                 preproc_y = standardize())

fit_plsc$singvals
fit_plsc$explained_cov

## ----scores-plot, fig.align='center'------------------------------------------
scores_df <- data.frame(
  LV1_x = scores(fit_plsc, "X")[, 1],
  LV1_y = scores(fit_plsc, "Y")[, 1],
  LV2_x = scores(fit_plsc, "X")[, 2],
  LV2_y = scores(fit_plsc, "Y")[, 2]
)

ggplot(scores_df, aes(x = LV1_y, y = LV1_x)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(x = "Behavior scores (LV1)", y = "Brain scores (LV1)",
       title = "Score association for LV1") +
  theme_minimal()

## ----perm-test----------------------------------------------------------------
set.seed(123)
pt <- perm_test(fit_plsc, X, Y, nperm = 199, comps = 3, parallel = FALSE)
pt$component_results
cat("Sequential n_significant (alpha = 0.05):", pt$n_significant, "\n")

## ----bootstrap----------------------------------------------------------------
set.seed(321)
boot_plsc <- bootstrap(fit_plsc, nboot = 120, X = X, Y = Y, comps = 2, parallel = FALSE)

# X-loadings bootstrap ratios for LV1
df_bsr <- data.frame(
  variable = paste0("x", seq_len(pX)),
  bsr      = boot_plsc$z_vx[, 1]
)

ggplot(df_bsr, aes(x = reorder(variable, bsr), y = bsr)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_col(fill = "#1f78b4") +
  coord_flip() +
  labs(x = NULL, y = "Bootstrap ratio (X loadings, LV1)",
       title = "Stable variables exceed |BSR| ≈ 2") +
  theme_minimal()

## ----bootstrap-y, echo=FALSE--------------------------------------------------
df_bsr_y <- data.frame(
  variable = paste0("y", seq_len(pY)),
  bsr      = boot_plsc$z_vy[, 1]
)

ggplot(df_bsr_y, aes(x = reorder(variable, bsr), y = bsr)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_col(fill = "#33a02c") +
  coord_flip() +
  labs(x = NULL, y = "Bootstrap ratio (Y loadings, LV1)",
       title = "Behavior variables driving LV1") +
  theme_minimal()

