# label switch: check all permutations
all_permutations <- function(x) {
  if (length(x) == 1L) {
    return(list(x))
  }
  out <- list()

  for (i in seq_along(x)) {
    rest <- x[-i]
    rest_perm <- all_permutations(rest)

    out <- c(
      out,
      lapply(rest_perm, function(p) c(x[i], p))
    )
  }

  out
}

# regression paths (to check the recovery of regression parameters)
get_regression_paths <- function(step2model) {
  parameter_table <- lavaan::lavaanify(step2model)

  regression_paths <- parameter_table[parameter_table$op == "~", c("lhs", "rhs")]

  regression_paths$path <- paste(
    regression_paths$lhs, "~", regression_paths$rhs)

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

  names(values) <- paths$path

  values
}

# step 2 evaluation
evaluate_step2_validation <- function(posteriors,
                                      beta_ks,
                                      true_cluster,
                                      true_beta,
                                      step2model,
                                      nclus = ncol(posteriors)) {

  # Convert posterior probabilities to hard cluster labels
  estimated_cluster <- max.col(posteriors, ties.method = "first")

  # Find best label mapping by trying all cluster-label permutations
  perms <- all_permutations(seq_len(nclus))

  misclass <- numeric(length(perms))
  relabeled_list <- vector("list", length(perms))

  for (i in seq_along(perms)) {
    # If perms[[i]] = c(2, 1), then:
    # estimated cluster 1 is treated as true cluster 2
    # estimated cluster 2 is treated as true cluster 1
    relabeled_list[[i]] <- perms[[i]][estimated_cluster]

    misclass[i] <- mean(relabeled_list[[i]] != true_cluster)
  }

  best_i <- which.min(misclass)

  best_mapping <- stats::setNames(
    perms[[best_i]],
    seq_len(nclus)
  )

  relabeled_cluster <- relabeled_list[[best_i]]

  # Extract model-specified regression paths
  paths <- get_regression_paths(step2model)

  # Compare beta estimates after cluster-label alignment
  beta_tables <- list()

  for (est_k in seq_len(nclus)) {
    true_k <- as.integer(unname(best_mapping[as.character(est_k)]))

    estimated_values <- extract_beta_values(
      beta_matrix = beta_ks[[est_k]],
      paths = paths
    )

    true_values <- extract_beta_values(
      beta_matrix = true_beta[, , true_k],
      paths = paths
    )

    difference <- estimated_values - true_values

    beta_tables[[est_k]] <- data.frame(
      estimated_cluster = est_k,
      matched_true_cluster = true_k,
      path = names(estimated_values),
      estimated = as.numeric(estimated_values),
      truth = as.numeric(true_values),
      difference = as.numeric(difference),
      stringsAsFactors = FALSE
    )
  }

  beta_table <- do.call(rbind, beta_tables)
  rownames(beta_table) <- NULL

  # Return useful summaries
  list(
    clustering = list(
      estimated_cluster = estimated_cluster,
      relabeled_cluster = relabeled_cluster,
      true_cluster = true_cluster,
      mapping = best_mapping,
      misclassification_error = misclass[best_i],
      correct_clustering = misclass[best_i] == 0,
      mean_uncertainty = mean(1 - apply(posteriors, 1, max))
    ),

    beta = list(
      paths = paths,
      table = beta_table,
      bias = mean(beta_table$difference),
      rmse = sqrt(mean(beta_table$difference^2)),
      max_abs_error = max(abs(beta_table$difference))
    )
  )
}


