get_signature_genes <- function(signature_id) {
  registry <- load_signature_registry()
  signatures <- attr(registry, "TBsignatures")
  split <- attr(registry, "TBsignaturesSplit")
  if (is.null(signatures) || !signature_id %in% names(signatures)) {
    return(list(genes = character(), up_genes = character(), down_genes = character(), has_direction = FALSE))
  }
  split_entry <- split[[signature_id]]
  if (is.list(split_entry)) {
    up_name <- grep("_up$", names(split_entry), value = TRUE)
    dn_name <- grep("_dn$|_down$", names(split_entry), value = TRUE)
    up <- unique(unlist(split_entry[up_name], use.names = FALSE))
    down <- unique(unlist(split_entry[dn_name], use.names = FALSE))
    return(list(genes = unique(c(up, down)), up_genes = up, down_genes = down, has_direction = TRUE))
  }
  list(genes = unique(signatures[[signature_id]]), up_genes = character(), down_genes = character(), has_direction = FALSE)
}

score_signature_signed_mean <- function(X, up_genes, down_genes) {
  up <- intersect(up_genes, rownames(X))
  down <- intersect(down_genes, rownames(X))
  up_score <- if (length(up)) colMeans(X[up, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(X))
  down_score <- if (length(down)) colMeans(X[down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(X))
  up_score - down_score
}

orient_signature_score <- function(score, y_dev, original_direction = NULL) {
  if (!is.null(original_direction) && original_direction == "higher_ATB") return(list(score = score, orientation = 1))
  auc <- compute_roc_auc(y_dev, score)
  if (!is.na(auc) && auc < 0.5) list(score = -score, orientation = -1) else list(score = score, orientation = 1)
}
