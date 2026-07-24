# Step 1 evaluation
# get the factor loadings paths
get_loading_paths <- function(step1model, obs_vars, markers = NULL) {
  par_table <- lavaan::lavaanify(step1model)

  loading_paths <- par_table[par_table$op == "=~", c("lhs", "rhs")]

  # Marker loadings are fixed to 1 during rescaling and are therefore excluded from the evaluation of estimated loadings
  if (!is.null(markers)) {
    marker_vars <- obs_vars[rowSums(markers == 1) > 0]
    loading_paths <- loading_paths[!loading_paths$rhs %in% marker_vars, , drop = FALSE]
  }

  rownames(loading_paths) <- NULL

  loading_paths
}

# extract factor loadings
extract_loading_values <- function(lambda_matrix, paths, obs_vars, lat_vars) {
  lambda_matrix <- as.matrix(lambda_matrix)

  # loadings are selected using observed and latent variable names
  rownames(lambda_matrix) <- obs_vars
  colnames(lambda_matrix) <- lat_vars

  values <- numeric(nrow(paths))

  for (i in seq_len(nrow(paths))) {
    latent <- paths$lhs[i]
    indicator <- paths$rhs[i]

    values[i] <- lambda_matrix[indicator, latent]
  }

  values
}

# Step 1 evaluation: recovery of MM-clusters, factor loadings, and factor covariance matrices
evaluate_step1_validation <- function(posteriors, lambda_ks, cov_eta, true_cluster,
                                      true_lambda, true_cov_eta, step1model,
                                      obs_vars, lat_vars, markers = NULL) {
  nclus <- ncol(posteriors)

  # 1. MM-cluster recovery
  # Cluster labels are arbitrary, so each possible relabeling is compared with the true cluster memberships
  estimated_cluster <- apply(posteriors, 1, which.max)
  permutations <- all_permutations(seq_len(nclus))
  misclass <- numeric(length(permutations))

  for (i in seq_along(permutations)) {
    relabeled_cluster <- permutations[[i]][estimated_cluster]
    misclass[i] <- mean(relabeled_cluster != true_cluster)
  }

  best_i <- which.min(misclass)
  best_mapping <- permutations[[best_i]]

  # 2. Loading recovery
  # The estimated clusters are first matched with the corresponding true clusters before their loadings are compared
  paths <- get_loading_paths(step1model = step1model, obs_vars = obs_vars, markers = markers)

  loading_difference <- numeric(0)

  for (est_k in seq_len(nclus)) {
    true_k <- best_mapping[est_k]

    estimated_values <- extract_loading_values(
      lambda_matrix = lambda_ks[[est_k]],
      paths = paths,
      obs_vars = obs_vars,
      lat_vars = lat_vars
    )

    true_values <- extract_loading_values(
      lambda_matrix = true_lambda[, , true_k],
      paths = paths,
      obs_vars = obs_vars,
      lat_vars = lat_vars
    )

    loading_difference <- c(loading_difference, estimated_values - true_values)
  }

  # 3. Factor covariance recovery
  # Only the lower triangle is used because the covariance matrices are symmetric
  lower_index <- lower.tri(cov_eta[, , 1], diag = TRUE)
  cov_difference <- numeric(0)

  for (g in seq_len(dim(cov_eta)[3])) {
    estimated_values <- cov_eta[, , g][lower_index]
    true_values <- true_cov_eta[, , g][lower_index]
    cov_difference <- c(cov_difference, estimated_values - true_values)
  }

  list(
    clustering = list(
      misclassification_error = misclass[best_i],
      mean_uncertainty = mean(1 - apply(posteriors, 1, max))),

    loading = list(
      rmse = sqrt(mean(loading_difference^2)),
      max_abs_error = max(abs(loading_difference))),

    factor_covariance = list(
      rmse = sqrt(mean(cov_difference^2)),
      max_abs_error = max(abs(cov_difference)))
  )
}
