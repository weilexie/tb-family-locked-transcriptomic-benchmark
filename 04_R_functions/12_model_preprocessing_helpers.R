# ML final classifier availability and matrix preparation
# Split from modeling_algorithm_combinations_final.R for maintainability.

automl_final_classifier_requirement <- function(classifier) {
  switch(classifier,
    ridge_logistic = "glmnet", lasso_logistic = "glmnet", elastic_net_025 = "glmnet", elastic_net_050 = "glmnet", elastic_net_075 = "glmnet",
    gbm = "gbm", xgboost = "xgboost", random_forest = "random_forest", lda = "MASS", naive_bayes = "e1071", pls_glm = "pls", svm_linear = "e1071", svm_rbf = "e1071", stepglm_logistic = "MASS", ""
  )
}

automl_final_classifier_available <- function(classifier) {
  req <- automl_final_classifier_requirement(classifier)
  if (identical(req, "")) return(TRUE)
  if (identical(req, "random_forest")) return(requireNamespace("ranger", quietly = TRUE) || requireNamespace("randomForest", quietly = TRUE))
  requireNamespace(req, quietly = TRUE)
}

automl_final_max_features_for_classifier <- function(classifier) {
  switch(classifier,
    svm_rbf = 80, svm_linear = 120, lda = 80, naive_bayes = 120, stepglm_logistic = 35, pls_glm = 120, gbm = 120, xgboost = 150, random_forest = 120, 250
  )
}

automl_final_prepare_xy <- function(X_train_raw, X_test_raw, train_ids, test_ids, selected, normalization, sig_train_all, sig_test_all, classifier) {
  gene_train <- matrix(nrow = length(train_ids), ncol = 0, dimnames = list(train_ids, character()))
  gene_test <- matrix(nrow = length(test_ids), ncol = 0, dimnames = list(test_ids, character()))
  if (length(selected$genes) > 0) {
    genes <- selected$genes[selected$genes %in% rownames(X_train_raw)]
    genes <- genes[rowMeans(!is.na(X_train_raw[genes, , drop = FALSE])) >= 0.50]
    if (length(genes) > 0) {
      norm_fit <- fit_normalizer(X_train_raw[genes, , drop = FALSE], normalization, genes)
      gene_train <- t(transform_normalizer(norm_fit, X_train_raw[genes, , drop = FALSE]))
      gene_test <- t(transform_normalizer(norm_fit, X_test_raw[genes, , drop = FALSE]))
      colnames(gene_train) <- paste0("gene__", make.names(colnames(gene_train), unique = TRUE))
      colnames(gene_test) <- colnames(gene_train)
    }
  }
  sig_train <- matrix(nrow = length(train_ids), ncol = 0, dimnames = list(train_ids, character()))
  sig_test <- matrix(nrow = length(test_ids), ncol = 0, dimnames = list(test_ids, character()))
  if (length(selected$signatures) > 0 && ncol(sig_train_all) > 0) {
    sigs <- selected$signatures[selected$signatures %in% colnames(sig_train_all)]
    if (length(sigs) > 0) {
      sig_train <- sig_train_all[train_ids, sigs, drop = FALSE]
      sig_test <- sig_test_all[test_ids, sigs, drop = FALSE]
      scaled <- automl_final_impute_scale_samples(sig_train, sig_test)
      sig_train <- scaled$train
      sig_test <- scaled$test
    }
  }
  x_train <- cbind(gene_train, sig_train)
  x_test <- cbind(gene_test, sig_test)
  max_features <- automl_final_max_features_for_classifier(classifier)
  reduced <- automl_final_reduce_zero_variance(x_train, x_test, max_features = max_features)
  reduced
}

