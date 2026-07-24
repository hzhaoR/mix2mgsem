# Get all possible cluster label orders
# This is needed because estimated cluster labels may be switched
all_permutations <- function(x) {
  if (length(x) == 1L) {return(list(x))}
  permutations <- list()

  for (i in seq_along(x)) {
    rest_perm <- all_permutations(x[-i]) # find all permutations of x with the current element removed

    for (j in seq_along(rest_perm)) {
      permutations[[length(permutations) + 1L]] <- c(x[i], rest_perm[[j]])
    }
  }

  permutations
}

# get the regression paths
get_regression_paths <- function(step2model) {
  par_table <- lavaan::lavaanify(step2model)

  regression_paths <- par_table[par_table$op == "~", c("lhs", "rhs")]
  rownames(regression_paths) <- NULL

  regression_paths
}

# extract betas from the matrix
extract_beta_values <- function(beta_matrix, paths) {
  values <- numeric(nrow(paths))

  for (i in seq_len(nrow(paths))) {
    lhs <- paths$lhs[i]
    rhs <- paths$rhs[i]

    values[i] <- beta_matrix[lhs, rhs]
  }

  values
}

# Step 2 evaluation: recovery of SR-clusters and regression parameters
evaluate_step2_validation <- function(posteriors, beta_ks, true_cluster,
                                      true_beta, step2model) {
  nclus <- ncol(posteriors)

  # 1. SR-cluster recovery
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

  # 2. Regression parameter recovery
  # The estimated clusters are matched with the corresponding true clusters before their beta matrices are compared
  paths <- get_regression_paths(step2model)

  beta_difference <- numeric(0)

  for (est_k in seq_len(nclus)) {
    true_k <- best_mapping[est_k]

    estimated_values <- extract_beta_values(
      beta_matrix = beta_ks[[est_k]],
      paths = paths
    )

    true_values <- extract_beta_values(
      beta_matrix = true_beta[, , true_k],
      paths = paths
    )

    beta_difference <- c(beta_difference, estimated_values - true_values)
  }

  list(
    clustering = list(
      misclassification_error = misclass[best_i],
      mean_uncertainty = mean(1 - apply(posteriors, 1, max))),

    beta = list(
      rmse = sqrt(mean(beta_difference^2)),
      max_abs_error = max(abs(beta_difference)))
  )
}


