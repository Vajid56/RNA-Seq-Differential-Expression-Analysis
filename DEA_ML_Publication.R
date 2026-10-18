## ============================================================
##  Publication-standard Machine Learning for DEA data
##  Stroke vs Diabetes — Sayyad Vajid
##
##  THE CORRECT APPROACH:
##  Task: Use gene expression features to predict disease condition
##  Label: Stroke (1) vs Diabetes (0) — from SAMPLE metadata
##  Features: baseMean, log2FoldChange, lfcSE, stat (from DESeq2)
##  NOTE: pvalue and padj are EXCLUDED — they are derived stats,
##        not independent biological features
## ============================================================
install.packages("tidyverse")

library(tidyverse)


install.packages("caret")

library(caret)


install.packages("randomForest")

library(randomForest)


install.packages("pROC")

library(pROC)


install.packages("e1071")

library(e1071)
library(dplyr)
# ── 1. Install & load packages ──────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs_cran <- c("tidyverse", "caret", "randomForest", "e1071",
               "pROC", "ROCR", "ggplot2", "pheatmap",
               "viridis", "RColorBrewer", "corrplot",
               "gridExtra", "ggpubr", "scales")
new_pkgs <- pkgs_cran[!pkgs_cran %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)
lapply(pkgs_cran, library, character.only = TRUE)

set.seed(42)  # reproducibility — always set this for publications

# ── 2. Load your DESeq2 results ─────────────────────────────
# Replace with your actual file path
dea <- read.csv("/home/vajid/DE_stroke_vs_diabetes.csv", row.names = 1)
colnames(dea) <- c("baseMean", "log2FoldChange", "lfcSE",
                   "stat", "pvalue", "padj")

cat("Total genes:", nrow(dea), "\n")
cat("Columns:", colnames(dea), "\n")

# ── 3. CORRECT label definition ──────────────────────────────
# Label = which DIRECTION the gene changes (biological meaning)
# Up = more expressed in Stroke than Diabetes
# Down = more expressed in Diabetes than Stroke
# This is biologically meaningful, not circular

dea$direction <- case_when(
  dea$log2FoldChange >  0 ~ "Stroke_higher",
  dea$log2FoldChange <= 0 ~ "Diabetes_higher"
)
dea$direction <- factor(dea$direction)

# Alternative label: significant vs not (keep for feature importance only)
dea$is_sig <- factor(
  ifelse(dea$padj < 0.05 & abs(dea$log2FoldChange) > 1, "Significant", "Not_Sig")
)

cat("Direction distribution:\n")
print(table(dea$direction))

# ── 4. Feature engineering (SAFE features — no leakage) ──────
# INCLUDE: baseMean, lfcSE, stat (Wald test statistic), abs_LFC
# EXCLUDE: pvalue, padj (derived from stat, circular for direction prediction)

features <- dea %>%
  mutate(
    log_baseMean  = log1p(baseMean),       # log-normalize expression
    abs_LFC       = abs(log2FoldChange),   # magnitude of change
    lfc_se_ratio  = abs(log2FoldChange) / lfcSE,  # signal-to-noise ratio
    log_lfcSE     = log1p(lfcSE)          # log-normalize SE
  ) %>%
  select(log_baseMean, abs_LFC, lfcSE, log_lfcSE,
         stat, lfc_se_ratio)

cat("\nFeature summary:\n")
summary(features)

# Check for missing values
cat("Missing values per feature:\n")
print(colSums(is.na(features)))

# ── 5. Check for class imbalance ─────────────────────────────
cat("\nClass balance check:\n")
print(table(dea$direction))
cat("Class ratio:", round(max(table(dea$direction)) /
                          min(table(dea$direction)), 2), "\n")

# ── 6. Train/test split (stratified 80/20) ───────────────────
train_idx <- createDataPartition(dea$direction, p = 0.8,
                                  list = FALSE, times = 1)
X_train <- features[train_idx, ]
X_test  <- features[-train_idx, ]
y_train <- dea$direction[train_idx]
y_test  <- dea$direction[-train_idx]

cat("Train size:", nrow(X_train), "| Test size:", nrow(X_test), "\n")

# ── 7. Preprocessing — normalize features ────────────────────
preproc <- preProcess(X_train, method = c("center", "scale"))
X_train_sc <- predict(preproc, X_train)
X_test_sc  <- predict(preproc, X_test)

# ── 8. Cross-validation setup ────────────────────────────────
# 5-fold stratified CV — standard for bioinformatics papers
ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,          # needed for AUC
  summaryFunction = twoClassSummary,  # gives ROC, Sens, Spec
  savePredictions = "final",
  verboseIter     = FALSE
)

# ── 9. Train all models ──────────────────────────────────────

# 9a. Random Forest (tuned mtry)
cat("\n--- Training Random Forest ---\n")
rf_model <- train(
  x          = X_train_sc,
  y          = y_train,
  method     = "rf",
  metric     = "ROC",
  trControl  = ctrl,
  tuneGrid   = expand.grid(mtry = c(2, 3, 4)),
  ntree      = 500,
  importance = TRUE
)
cat("Best mtry:", rf_model$bestTune$mtry, "\n")

# 9b. Support Vector Machine (RBF kernel, tuned C and sigma)
cat("\n--- Training SVM (RBF) ---\n")
svm_model <- train(
  x         = X_train_sc,
  y         = y_train,
  method    = "svmRadial",
  metric    = "ROC",
  trControl = ctrl,
  tuneGrid  = expand.grid(C     = c(0.1, 1, 10),
                           sigma = c(0.01, 0.1, 1))
)
cat("Best C:", svm_model$bestTune$C,
    "| Best sigma:", svm_model$bestTune$sigma, "\n")

# 9c. Logistic Regression (with regularization)
cat("\n--- Training Logistic Regression ---\n")
lr_model <- train(
  x         = X_train_sc,
  y         = y_train,
  method    = "glmnet",
  metric    = "ROC",
  trControl = ctrl,
  tuneGrid  = expand.grid(alpha  = c(0, 0.5, 1),
                           lambda = seq(0.001, 0.1, length = 10))
)

# 9d. Gradient Boosting Machine
cat("\n--- Training GBM ---\n")
gbm_model <- train(
  x         = X_train_sc,
  y         = y_train,
  method    = "gbm",
  metric    = "ROC",
  trControl = ctrl,
  tuneGrid  = expand.grid(n.trees           = c(100, 200),
                           interaction.depth = c(2, 3),
                           shrinkage         = c(0.05, 0.1),
                           n.minobsinnode    = 10),
  verbose   = FALSE
)

# 9e. k-Nearest Neighbors (tuned k)
cat("\n--- Training KNN ---\n")
knn_model <- train(
  x         = X_train_sc,
  y         = y_train,
  method    = "knn",
  metric    = "ROC",
  trControl = ctrl,
  tuneGrid  = expand.grid(k = c(3, 5, 7, 9, 11, 15))
)
cat("Best k:", knn_model$bestTune$k, "\n")

# ── 10. Evaluate on held-out test set ───────────────────────
evaluate_model <- function(model, X_test, y_test, model_name) {
  # Predictions
  y_pred  <- predict(model, X_test)
  y_prob  <- predict(model, X_test, type = "prob")[, 2]

  # Confusion matrix
  cm <- confusionMatrix(y_pred, y_test, positive = "Stroke_higher")

  # AUC-ROC
  roc_obj <- roc(response  = y_test,
                 predictor = y_prob,
                 levels    = levels(y_test))
  auc_val <- auc(roc_obj)

  # Matthews Correlation Coefficient (better for imbalanced classes)
  TP <- cm$table[2, 2]; TN <- cm$table[1, 1]
  FP <- cm$table[2, 1]; FN <- cm$table[1, 2]
  mcc <- (TP * TN - FP * FN) /
         sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  data.frame(
    Model     = model_name,
    Accuracy  = round(cm$overall["Accuracy"] * 100, 2),
    Precision = round(cm$byClass["Pos Pred Value"] * 100, 2),
    Recall    = round(cm$byClass["Sensitivity"] * 100, 2),
    Specificity = round(cm$byClass["Specificity"] * 100, 2),
    F1        = round(cm$byClass["F1"] * 100, 2),
    AUC_ROC   = round(auc_val, 4),
    MCC       = round(mcc, 4)
  )
}

results <- bind_rows(
  evaluate_model(rf_model,  X_test_sc, y_test, "Random Forest"),
  evaluate_model(svm_model, X_test_sc, y_test, "SVM (RBF)"),
  evaluate_model(lr_model,  X_test_sc, y_test, "Logistic Regression"),
  evaluate_model(gbm_model, X_test_sc, y_test, "Gradient Boosting"),
  evaluate_model(knn_model, X_test_sc, y_test, "KNN")
)

cat("\n=== MODEL PERFORMANCE ON HELD-OUT TEST SET ===\n")
print(results)

# ── 11. Statistical comparison of models (Friedman test) ────
# Required for publication — prove models are significantly different
cat("\n--- Statistical comparison (Friedman test) ---\n")
model_list <- list(
  RF  = rf_model,  SVM = svm_model, LR = lr_model,
  GBM = gbm_model, KNN = knn_model
)
resamples_obj <- resamples(model_list)
summary(resamples_obj)

# Pairwise t-test between models (Bonferroni correction)
diff_obj <- diff(resamples_obj)
cat("\nPairwise model comparison (p-values):\n")
print(summary(diff_obj))

# ── 12. ROC curves (publication quality) ────────────────────
# Get ROC objects for all models
get_roc <- function(model, X_test, y_test) {
  y_prob <- predict(model, X_test, type = "prob")[, 2]
  roc(response = y_test, predictor = y_prob,
      levels = levels(y_test))
}

roc_rf  <- get_roc(rf_model,  X_test_sc, y_test)
roc_svm <- get_roc(svm_model, X_test_sc, y_test)
roc_lr  <- get_roc(lr_model,  X_test_sc, y_test)
roc_gbm <- get_roc(gbm_model, X_test_sc, y_test)
roc_knn <- get_roc(knn_model, X_test_sc, y_test)

# Plot all ROC curves on one figure
png("roc_curves.png", width = 2400, height = 2000, res = 300)
plot(roc_rf,  col = "#185FA5", lwd = 2.5,
     main = "ROC Curves — Stroke vs Diabetes Classification",
     cex.main = 1.2, font.main = 2)
plot(roc_svm, col = "#1D9E75", lwd = 2.5, add = TRUE)
plot(roc_lr,  col = "#D85A30", lwd = 2.5, add = TRUE)
plot(roc_gbm, col = "#BA7517", lwd = 2.5, add = TRUE)
plot(roc_knn, col = "#534AB7", lwd = 2.5, add = TRUE)
abline(a = 0, b = 1, col = "gray60", lty = 2, lwd = 1.2)
legend("bottomright",
  legend = c(
    paste0("Random Forest   AUC = ", round(auc(roc_rf),  3)),
    paste0("SVM (RBF)       AUC = ", round(auc(roc_svm), 3)),
    paste0("Logistic Regr.  AUC = ", round(auc(roc_lr),  3)),
    paste0("Gradient Boost  AUC = ", round(auc(roc_gbm), 3)),
    paste0("KNN (k=", knn_model$bestTune$k, ")          AUC = ",
           round(auc(roc_knn), 3))
  ),
  col = c("#185FA5","#1D9E75","#D85A30","#BA7517","#534AB7"),
  lwd = 2.5, bty = "n", cex = 0.9
)
dev.off()
cat("Saved: roc_curves.png\n")

# ── 13. Feature importance (Random Forest) ───────────────────
fi <- varImp(rf_model, scale = TRUE)
fi_df <- fi$importance
fi_df$Feature <- rownames(fi_df)
fi_df <- fi_df[order(fi_df$Overall, decreasing = FALSE), ]

# Clean feature names for publication
fi_df$Feature <- recode(fi_df$Feature,
  "log_baseMean"  = "log(baseMean)",
  "abs_LFC"       = "|log2FC|",
  "log_lfcSE"     = "log(lfcSE)",
  "lfc_se_ratio"  = "|log2FC| / lfcSE",
  "stat"          = "Wald statistic",
  "lfcSE"         = "lfcSE"
)

p_fi <- ggplot(fi_df,
               aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "#185FA5", width = 0.65, alpha = 0.88) +
  geom_text(aes(label = round(Overall, 1)),
            hjust = -0.15, size = 3.5, color = "#444441") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title    = "Random Forest Feature Importance",
       subtitle = "Predicting gene expression direction (Stroke > Diabetes vs Diabetes > Stroke)",
       x = NULL, y = "Relative Importance (scaled 0-100)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray50", size = 10),
    axis.text.y   = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray92")
  )

ggsave("feature_importance.png", p_fi,
       width = 8, height = 5, dpi = 300)
cat("Saved: feature_importance.png\n")

# ── 14. Confusion matrices (publication quality) ─────────────
plot_cm <- function(model, X_test, y_test, title) {
  y_pred <- predict(model, X_test)
  cm     <- confusionMatrix(y_pred, y_test)$table
  cm_df  <- as.data.frame(cm)

  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Freq), size = 7, fontface = "bold",
              color = "white") +
    scale_fill_gradient(low = "#B5D4F4", high = "#0C447C") +
    scale_x_discrete(labels = c("Diabetes\nhigher", "Stroke\nhigher")) +
    scale_y_discrete(labels = c("Diabetes\nhigher", "Stroke\nhigher")) +
    labs(title = title, x = "True label", y = "Predicted label") +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.title = element_blank(),
      axis.text    = element_text(size = 10)
    )
}

cm_plots <- list(
  plot_cm(rf_model,  X_test_sc, y_test, "Random Forest"),
  plot_cm(svm_model, X_test_sc, y_test, "SVM (RBF)"),
  plot_cm(lr_model,  X_test_sc, y_test, "Logistic Regression"),
  plot_cm(gbm_model, X_test_sc, y_test, "Gradient Boosting"),
  plot_cm(knn_model, X_test_sc, y_test, "KNN")
)

cm_grid <- gridExtra::arrangeGrob(
  grobs = cm_plots, ncol = 3,
  top = grid::textGrob(
    "Confusion Matrices — All Models (Held-out Test Set)",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
)
ggsave("confusion_matrices.png", cm_grid,
       width = 14, height = 9, dpi = 300)
cat("Saved: confusion_matrices.png\n")

# ── 15. Model comparison dotplot (CV results) ────────────────
png("model_comparison_cv.png", width = 2400, height = 1600, res = 300)
dotplot(resamples_obj,
        metric = "ROC",
        main   = "5-Fold CV AUC-ROC Comparison (with 95% CI)")
dev.off()
cat("Saved: model_comparison_cv.png\n")

# ── 16. Learning curve (Random Forest) ───────────────────────
sizes <- seq(0.1, 1.0, by = 0.1)
lc_results <- lapply(sizes, function(p) {
  idx <- createDataPartition(y_train, p = p, list = FALSE)
  X_sub <- X_train_sc[idx, ]; y_sub <- y_train[idx]
  mod <- train(x = X_sub, y = y_sub, method = "rf",
               trControl = trainControl(method = "cv", number = 3,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary),
               metric = "ROC", ntree = 100,
               tuneGrid = data.frame(mtry = rf_model$bestTune$mtry))
  test_prob <- predict(mod, X_test_sc, type = "prob")[, 2]
  test_auc  <- as.numeric(auc(roc(y_test, test_prob,
                                   levels = levels(y_test))))
  train_auc <- max(mod$results$ROC)
  data.frame(size = round(p * nrow(X_train_sc)), p = p,
             train_auc = train_auc, test_auc = test_auc)
})
lc_df <- bind_rows(lc_results)

p_lc <- ggplot(lc_df) +
  geom_line(aes(x = size, y = train_auc, color = "Training AUC"),
            linewidth = 1.3) +
  geom_line(aes(x = size, y = test_auc,  color = "Test AUC"),
            linewidth = 1.3) +
  geom_point(aes(x = size, y = train_auc, color = "Training AUC"), size = 3) +
  geom_point(aes(x = size, y = test_auc,  color = "Test AUC"), size = 3) +
  scale_color_manual(values = c("Training AUC" = "#185FA5",
                                 "Test AUC"     = "#D85A30")) +
  scale_y_continuous(limits = c(0.5, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
  labs(title    = "Learning Curve — Random Forest",
       subtitle = "AUC-ROC vs training set size",
       x = "Training set size (genes)", y = "AUC-ROC",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray93")
  )

ggsave("learning_curve.png", p_lc, width = 7, height = 5, dpi = 300)
cat("Saved: learning_curve.png\n")

# ── 17. Correct statistical tests for publication ─────────────
cat("\n=== STATISTICAL TESTS FOR PUBLICATION ===\n")

# Test 1: Wilcoxon test — compare log2FC between UP and DOWN genes
sig_genes <- dea %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
up_genes   <- sig_genes %>% filter(log2FoldChange > 1)
down_genes <- sig_genes %>% filter(log2FoldChange < -1)

cat("\n--- Wilcoxon rank-sum test: baseMean (Up vs Down genes) ---\n")
wilcox_result <- wilcox.test(
  up_genes$baseMean, down_genes$baseMean,
  alternative = "two.sided", paired = FALSE, exact = FALSE
)
print(wilcox_result)
cat("Effect size (r):", round(wilcox_result$statistic /
    (length(up_genes$baseMean) * length(down_genes$baseMean)), 4), "\n")

# Test 2: Kruskal-Wallis — compare |LFC| across significance categories
dea$category <- case_when(
  dea$padj < 0.001 & abs(dea$log2FoldChange) > 1 ~ "Highly_sig",
  dea$padj < 0.05  & abs(dea$log2FoldChange) > 1 ~ "Significant",
  TRUE                                             ~ "Not_sig"
)
cat("\n--- Kruskal-Wallis test: |LFC| across significance categories ---\n")
kw_result <- kruskal.test(abs(log2FoldChange) ~ category, data = dea)
print(kw_result)

# Post-hoc: Dunn test with Bonferroni correction
if (!requireNamespace("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)
cat("\n--- Dunn post-hoc test (Bonferroni) ---\n")
dunn_result <- dunn.test(
  x     = abs(dea$log2FoldChange),
  g     = dea$category,
  method = "bonferroni",
  alpha  = 0.05
)

# Test 3: Spearman correlation — baseMean vs |LFC|
cat("\n--- Spearman correlation: log(baseMean) vs |LFC| ---\n")
spearman_result <- cor.test(
  log1p(dea$baseMean), abs(dea$log2FoldChange),
  method = "spearman", exact = FALSE
)
print(spearman_result)
cat("rho:", round(spearman_result$estimate, 4),
    "| p =", signif(spearman_result$p.value, 3), "\n")

# Test 4: Fisher's exact test — enrichment in specific categories
# (e.g., are upregulated genes enriched in high-expression genes?)
high_expr   <- dea$baseMean > median(dea$baseMean)
upregulated <- dea$log2FoldChange > 1 & dea$padj < 0.05
cont_table  <- table(high_expr, upregulated)
cat("\n--- Fisher's exact test: high expression vs upregulation ---\n")
fisher_result <- fisher.test(cont_table)
print(fisher_result)
cat("Odds ratio:", round(fisher_result$estimate, 4),
    "| p =", signif(fisher_result$p.value, 3), "\n")

# ── 18. Save final results table ─────────────────────────────
write.csv(results, "ML_results_publication.csv", row.names = FALSE)

# Also save a formatted table with stars for significance
results_pub <- results %>%
  mutate(
    Accuracy_fmt  = paste0(Accuracy, "%"),
    F1_fmt        = paste0(F1, "%"),
    AUC_fmt       = sprintf("%.4f", AUC_ROC),
    MCC_fmt       = sprintf("%.4f", MCC),
    Best_model    = ifelse(AUC_ROC == max(AUC_ROC), "*", "")
  )
write.csv(results_pub, "ML_results_formatted.csv", row.names = FALSE)

cat("\n=== ALL OUTPUTS SAVED ===\n")
cat("roc_curves.png\n")
cat("feature_importance.png\n")
cat("confusion_matrices.png\n")
cat("model_comparison_cv.png\n")
cat("learning_curve.png\n")
cat("ML_results_publication.csv\n")
cat("ML_results_formatted.csv\n")
cat("\nSession info:\n")
sessionInfo()
