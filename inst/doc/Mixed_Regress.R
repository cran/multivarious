params <-
list(family = "red")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
library(multivarious)

## ----simulate_data------------------------------------------------------------
set.seed(1)

n_subject <- 18
levels_within <- c("low", "mid", "high")

design <- expand.grid(
  subject = factor(seq_len(n_subject)),
  level = factor(levels_within, levels = levels_within),
  KEEP.OUT.ATTRS = FALSE
)

subject_group <- rep(c("A", "B"), length.out = n_subject)
design$group <- factor(subject_group[as.integer(design$subject)])

level_num <- c(low = -1, mid = 0, high = 1)[as.character(design$level)]
group_num <- ifelse(design$group == "B", 1, 0)
subj_idx <- as.integer(design$subject)

b0 <- rnorm(n_subject, sd = 0.7)
b1 <- rnorm(n_subject, sd = 0.3)

n <- nrow(design)
p <- 8

Y <- cbind(
  b0[subj_idx] + 1.2 * level_num + rnorm(n, sd = 0.2),
  0.8 * group_num + rnorm(n, sd = 0.2),
  1.4 * level_num * group_num + rnorm(n, sd = 0.2),
  -0.9 * level_num + rnorm(n, sd = 0.2),
  b1[subj_idx] * level_num + rnorm(n, sd = 0.2),
  rnorm(n, sd = 0.2),
  rnorm(n, sd = 0.2),
  rnorm(n, sd = 0.2)
)

dim(Y)

## ----fit_model----------------------------------------------------------------
fit <- mixed_regress(
  Y,
  design = design,
  fixed = ~ group * level,
  random = ~ 1 + level | subject,
  basis = shared_pca(ncomp = 4),
  preproc = center()
)

print(fit)
summary(fit)

## ----extract_effect-----------------------------------------------------------
E <- effect(fit, "group:level")

print(E)
ncomp(E)
components(E)[1:8, ]

## ----reconstruct_effect-------------------------------------------------------
E_proc <- reconstruct(E, scale = "processed")
E_orig <- reconstruct(E, scale = "original")

dim(E_proc)
dim(E_orig)
round(E_orig[1:6, 1:4], 3)

## ----permutation_test---------------------------------------------------------
set.seed(2)
pt <- perm_test(E, nperm = 99, alpha = 0.10)

print(pt)
pt$component_results

## ----truncate_to_selected_rank------------------------------------------------
k <- ncomp(pt)
E_sig <- truncate(E, k)

k
ncomp(E_sig)

## ----bootstrap_effect---------------------------------------------------------
set.seed(3)
bres <- bootstrap(E, nboot = 49, resample = "subject")

print(bres)
bres$singular_values_mean

## ----array_input--------------------------------------------------------------
Y_array <- array(NA_real_, dim = c(n_subject, length(levels_within), p))
idx <- 1
for (i in seq_len(n_subject)) {
  for (j in seq_along(levels_within)) {
    Y_array[i, j, ] <- Y[idx, ]
    idx <- idx + 1
  }
}

fit_array <- mixed_regress(
  Y_array,
  design = design,
  fixed = ~ group * level,
  random = ~ 1 | subject,
  basis = shared_pca(ncomp = 4),
  preproc = center()
)

effect(fit_array, "level")$term

## ----internal_checks, eval=nzchar(Sys.getenv("_MULTIVARIOUS_DEV_COVERAGE")), include=FALSE----
set.seed(4)
chk_fit <- mixed_regress(
  Y,
  design = design,
  fixed = ~ group * level,
  random = ~ 1 | subject,
  basis = shared_pca(ncomp = 3),
  preproc = pass()
)
chk_E <- effect(chk_fit, "group:level")
chk_pt <- perm_test(chk_E, nperm = 19, alpha = 0.2)
stopifnot(inherits(chk_E, "effect_operator"))
stopifnot(is.numeric(chk_pt$omnibus_p_value))

