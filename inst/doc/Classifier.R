params <-
list(family = "red")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4
)
library(multivarious)
library(dplyr)
library(tibble)
library(ggplot2) # Needed for plots
library(knitr)   # Needed for kable


## ----data_classifier----------------------------------------------------------
data(iris)
X   <- as.matrix(iris[, 1:4])
grp <- iris$Species

# Fit classical Linear DA and wrap it
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("MASS package required for LDA example")
}


# 1. Define and fit the pre-processing step using the training data
preproc_fitted <- fit(center(), X)
# 2. Transform the data
Xp <- transform(preproc_fitted, X)

# Assuming discriminant_projector, prep, center, scores are available
lda_fit <- MASS::lda(X, grouping = grp)


disc_proj <- multivarious::discriminant_projector(
  v      = lda_fit$scaling,                 # loadings (p × d)
  s      = Xp %*% lda_fit$scaling,           # scores   (n × d)
  sdev   = lda_fit$svd,                     # singular values
  labels = grp,
  preproc = preproc_fitted            # Pass the fitted pre-processor
)
print(disc_proj)

## ----plot_latent_space--------------------------------------------------------
scores_df <- as_tibble(scores(disc_proj)[, 1:2],
                       .name_repair = ~ c("LD1","LD2")) |>
  mutate(Species = iris$Species)

ggplot(scores_df, aes(LD1, LD2, colour = Species)) +
  geom_point(size = 2, alpha = .7) +
  stat_ellipse(level = .9, linewidth = .3) +
  theme_minimal() +
  ggtitle("Iris – first two LDA components")

## ----build_knn_classifier-----------------------------------------------------
set.seed(42)
train_id <- sample(seq_len(nrow(X)), size = 0.7*nrow(X))
test_id  <- setdiff(seq_len(nrow(X)), train_id)

# Assuming classifier function is available
clf_knn <- multivarious::classifier(
  x       = disc_proj,
  labels  = grp[train_id],
  new_data= X[train_id, ],  # Use training data to get reference scores
  knn     = 3
)
print(clf_knn)

## ----predict_knn--------------------------------------------------------------
pred_knn <- predict(clf_knn, new_data = X[test_id, ],
                    metric = "cosine", prob_type = "knn_proportion")

head(pred_knn$prob, 3)
print(paste("Overall Accuracy:", mean(pred_knn$class == grp[test_id])))

# Assuming rank_score and topk are available
rk  <- rank_score(pred_knn$prob, grp[test_id])
tk2 <- topk      (pred_knn$prob, grp[test_id], k = 2)

tibble(
  prank_mean = mean(rk$prank),
  top2_acc   = mean(tk2$topk)
)

## ----plot_confusion_matrix----------------------------------------------------
cm <- table(
  Truth     = grp[test_id],
  Predicted = pred_knn$class
)

# Heat-map
cm_df <- as.data.frame(cm)
ggplot(cm_df, aes(Truth, Predicted, fill = Freq)) +
  geom_tile(colour = "grey80") +
  geom_text(aes(label = Freq), colour = "white", size = 4) +
  scale_fill_gradient(low = "#4575b4", high = "#d73027", name="Count", limits = c(0, 15)) +
  scale_y_discrete(limits = rev(levels(cm_df$Predicted))) +
  theme_minimal(base_size = 12) + coord_equal() +
  ggtitle("k-NN (k = 3) confusion matrix – test set") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Pretty table as well
knitr::kable(cm, caption = "Confusion matrix (counts)")

## ----build_rf_classifier------------------------------------------------------
# Check if randomForest is installed
if (requireNamespace("randomForest", quietly = TRUE)) {

  # Assuming rf_classifier.projector method is available
  rf_clf <- rf_classifier( # Using the generic here
    x       = disc_proj,
    labels  = grp[train_id],
    # Pass scores directly if method requires it, or let it call scores(x)
    scores  = scores(disc_proj)[train_id, ]
  )

  pred_rf <- predict(rf_clf, new_data = X[test_id, ])
  print(paste("RF Accuracy:", mean(pred_rf$class == grp[test_id])))
} else {
  cat("randomForest package not installed. Skipping RF example.\n")
}


## ----predict_partial----------------------------------------------------------
sepal_cols <- 1:2

# Create a classifier using reference scores from Sepal columns only
clf_knn_sepal <- multivarious::classifier(
  x       = disc_proj,
  labels  = grp[train_id],
  new_data= X[train_id, sepal_cols],  # Use training data subset
  colind  = sepal_cols,             # Indicate which columns were used
  knn     = 3
)

# Predict using the dedicated sepal classifier
pred_sepal <- predict(
  clf_knn_sepal,                  # Use the sepal-specific classifier
  new_data = X[test_id, sepal_cols]
  # No need for colind here as clf_knn_sepal expects sepal data
)

print(paste("Accuracy (Sepal only):", mean(pred_sepal$class == grp[test_id])))

## ----calc_feature_importance--------------------------------------------------
blocks <- list(
  Sepal = 1:2,
  Petal = 3:4
)

# Assuming feature_importance is available
fi <- feature_importance(
  clf_knn,
  new_data  = X[test_id, ],
  true_labels = grp[test_id], # Pass the correct test set labels
  blocks    = blocks,
  fun       = rank_score, # Use rank_score as the performance metric
  fun_direction = "lower_is_better",
  approach  = "marginal" # Calculate marginal drop when block is removed
)
print(fi)

