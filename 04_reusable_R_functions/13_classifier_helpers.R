# ML final classifier fitting and prediction
# Split from modeling_algorithm_combinations_final.R for maintainability.

automl_final_fit_predict <- function(classifier, x_train, y_train_label, x_test, seed = 20260501) {
  if (ncol(x_train) < 1) stop("no_features", call. = FALSE)
  y <- as_binary(y_train_label)
  if (length(unique(y)) < 2) stop("training_set_has_one_class", call. = FALSE)
  set.seed(seed)
  if (classifier %in% c("ridge_logistic", "lasso_logistic", "elastic_net_025", "elastic_net_050", "elastic_net_075")) {
    alpha <- switch(classifier, ridge_logistic = 0, lasso_logistic = 1, elastic_net_025 = 0.25, elastic_net_050 = 0.50, elastic_net_075 = 0.75)
    fit <- glmnet::glmnet(x_train, y, family = "binomial", alpha = alpha, weights = class_weights(y))
    pred_train <- stats::predict(fit, newx = x_train, type = "response")
    aucs <- apply(pred_train, 2, function(p) compute_roc_auc(y, as.numeric(p)))
    dfs <- colSums(abs(as.matrix(stats::coef(fit))[-1, , drop = FALSE]) > 1e-8)
    max_auc <- max(aucs, na.rm = TRUE)
    max_df <- if (classifier == "ridge_logistic") Inf else 80
    eligible <- which(aucs >= max_auc - 0.01 & dfs <= max_df)
    if (length(eligible) == 0) eligible <- which.max(ifelse(is.na(aucs), -Inf, aucs))
    idx <- eligible[which.min(dfs[eligible])]
    return(as.numeric(stats::predict(fit, newx = x_test, s = fit$lambda[idx], type = "response")))
  }
  if (classifier == "stepglm_logistic") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    df$y <- y
    fit0 <- suppressWarnings(stats::glm(y ~ 1, data = df, family = stats::binomial()))
    fit1 <- suppressWarnings(stats::glm(y ~ ., data = df, family = stats::binomial()))
    st <- suppressWarnings(stats::step(fit0, scope = list(lower = fit0, upper = fit1), direction = "both", trace = 0, steps = 15))
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- make.names(names(tx), unique = TRUE)
    return(as.numeric(stats::predict(st, newdata = tx, type = "response")))
  }
  if (classifier == "lda") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    df$y <- factor(ifelse(y == 1, "ATB", "nonATB"), levels = c("nonATB", "ATB"))
    fit <- MASS::lda(y ~ ., data = df)
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- make.names(names(tx), unique = TRUE)
    return(as.numeric(stats::predict(fit, tx)$posterior[, "ATB"]))
  }
  if (classifier == "naive_bayes") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    yfac <- factor(ifelse(y == 1, "ATB", "nonATB"), levels = c("nonATB", "ATB"))
    fit <- e1071::naiveBayes(df, yfac)
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- names(df)
    pp <- stats::predict(fit, tx, type = "raw")
    return(as.numeric(pp[, "ATB"]))
  }
  if (classifier == "pls_glm") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    df$y <- y
    ncomp <- min(5, ncol(x_train), max(1, nrow(x_train) - 2))
    fit <- pls::plsr(y ~ ., data = df, ncomp = ncomp, validation = "none")
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- make.names(colnames(x_train), unique = TRUE)
    score <- as.numeric(stats::predict(fit, newdata = tx, ncomp = ncomp))
    return(pmin(pmax(score, 0), 1))
  }
  if (classifier %in% c("svm_linear", "svm_rbf")) {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    yfac <- factor(ifelse(y == 1, "ATB", "nonATB"), levels = c("nonATB", "ATB"))
    if (classifier == "svm_linear") {
      fit <- e1071::svm(df, yfac, kernel = "linear", probability = TRUE, scale = FALSE, cost = 1)
    } else {
      fit <- e1071::svm(df, yfac, kernel = "radial", probability = TRUE, scale = FALSE, cost = 1, gamma = 1 / max(ncol(x_train), 1))
    }
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- names(df)
    pred <- stats::predict(fit, tx, probability = TRUE)
    prob <- attr(pred, "probabilities")
    if (is.null(prob) || !"ATB" %in% colnames(prob)) return(as.numeric(pred == "ATB"))
    return(as.numeric(prob[, "ATB"]))
  }
  if (classifier == "gbm") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    df$y <- y
    fit <- suppressWarnings(gbm::gbm(y ~ ., data = df, distribution = "bernoulli", n.trees = 70, interaction.depth = 2, shrinkage = 0.05, bag.fraction = 0.8, train.fraction = 1, verbose = FALSE))
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- make.names(colnames(x_train), unique = TRUE)
    return(as.numeric(stats::predict(fit, newdata = tx, n.trees = 70, type = "response")))
  }
  if (classifier == "xgboost") {
    fit <- xgboost::xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2, eta = 0.06, subsample = 0.8, colsample_bytree = 0.8, nthread = 1), data = xgboost::xgb.DMatrix(x_train, label = y), nrounds = 70, verbose = 0)
    return(as.numeric(stats::predict(fit, xgboost::xgb.DMatrix(x_test))))
  }
  if (classifier == "random_forest") {
    df <- as.data.frame(x_train, check.names = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    df$y <- factor(ifelse(y == 1, "ATB", "nonATB"), levels = c("nonATB", "ATB"))
    if (requireNamespace("ranger", quietly = TRUE)) {
      fit <- ranger::ranger(y ~ ., data = df, probability = TRUE, num.trees = 120, mtry = max(1, floor(sqrt(ncol(x_train)))), seed = seed, num.threads = 1)
      tx <- as.data.frame(x_test, check.names = FALSE)
      names(tx) <- make.names(colnames(x_train), unique = TRUE)
      pp <- stats::predict(fit, tx)$predictions
      return(as.numeric(pp[, "ATB"]))
    }
    fit <- randomForest::randomForest(x = df[, setdiff(names(df), "y"), drop = FALSE], y = df$y, ntree = 120)
    tx <- as.data.frame(x_test, check.names = FALSE)
    names(tx) <- names(df)[names(df) != "y"]
    pp <- stats::predict(fit, tx, type = "prob")
    return(as.numeric(pp[, "ATB"]))
  }
  stop(sprintf("unknown_classifier_%s", classifier), call. = FALSE)
}

