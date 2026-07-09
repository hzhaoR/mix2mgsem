# Step 1 validation

# label switch: check all permutations
# all_permutations <- function(x) {
#   if (length(x) == 1L) {
#     return(list(x))
#   }
#
#   out <- list()
#
#   for (i in seq_along(x)) {
#     rest <- x[-i]
#     rest_perm <- all_permutations(rest)
#
#     out <- c(
#       out,
#       lapply(rest_perm, function(p) c(x[i], p))
#     )
#   }
#
#   out
# }

# factor loadings
get_loading_paths <- function(step1model,
                              markers = NULL,
                              include_marker_loadings = FALSE) {
  parameter_table <- lavaan::lavaanify(step1model)

  loading_paths <- parameter_table[
    parameter_table$op == "=~",
    c("lhs", "rhs")
  ]

  loading_paths$lhs <- as.character(loading_paths$lhs)
  loading_paths$rhs <- as.character(loading_paths$rhs)

  loading_paths <- unique(loading_paths)

  if (!include_marker_loadings && !is.null(markers)) {
    markers <- as.character(unlist(markers))
    loading_paths <- loading_paths[
      !loading_paths$rhs %in% markers,
      ,
      drop = FALSE
    ]
  }

  loading_paths$path <- paste(
    loading_paths$lhs,
    "=~",
    loading_paths$rhs
  )

  rownames(loading_paths) <- NULL

  loading_paths
}

# extract factor loadings
extract_loading_values <- function(lambda_matrix,
                                   paths,
                                   observed_vars,
                                   latent_vars) {
  lambda_matrix <- as.matrix(lambda_matrix)

  rownames(lambda_matrix) <- observed_vars
  colnames(lambda_matrix) <- latent_vars

  values <- numeric(nrow(paths))

  for (i in seq_len(nrow(paths))) {
    latent <- paths$lhs[i]
    indicator <- paths$rhs[i]

    values[i] <- lambda_matrix[indicator, latent]
  }

  names(values) <- paths$path

  values
}

# step 1 evaluation
evaluate_step1_validation <- function(posteriors,
                                      lambda_ks,
                                      cov_eta,
                                      true_cluster,
                                      true_lambda,
                                      true_cov_eta,
                                      step1model,
                                      observed_vars,
                                      latent_vars,
                                      markers = NULL,
                                      include_marker_loadings = FALSE,
                                      nclus = ncol(posteriors)) {

  # ---------------------------------------------------------------------------
  # 1. MM-cluster recovery
  # ---------------------------------------------------------------------------

  estimated_cluster <- max.col(posteriors, ties.method = "first")

  perms <- all_permutations(seq_len(nclus))

  misclass <- numeric(length(perms))
  relabeled_list <- vector("list", length(perms))

  for (i in seq_along(perms)) {
    relabeled_list[[i]] <- perms[[i]][estimated_cluster]
    misclass[i] <- mean(relabeled_list[[i]] != true_cluster)
  }

  best_i <- which.min(misclass)

  best_mapping <- stats::setNames(
    perms[[best_i]],
    seq_len(nclus)
  )

  relabeled_cluster <- relabeled_list[[best_i]]

  # ---------------------------------------------------------------------------
  # 2. Loading recovery
  # ---------------------------------------------------------------------------

  paths <- get_loading_paths(
    step1model = step1model,
    markers = markers,
    include_marker_loadings = include_marker_loadings
  )

  loading_tables <- list()

  for (est_k in seq_len(nclus)) {
    true_k <- as.integer(unname(best_mapping[as.character(est_k)]))

    estimated_values <- extract_loading_values(
      lambda_matrix = lambda_ks[[est_k]],
      paths = paths,
      observed_vars = observed_vars,
      latent_vars = latent_vars
    )

    true_values <- extract_loading_values(
      lambda_matrix = true_lambda[, , true_k],
      paths = paths,
      observed_vars = observed_vars,
      latent_vars = latent_vars
    )

    difference <- estimated_values - true_values

    loading_tables[[est_k]] <- data.frame(
      estimated_cluster = est_k,
      matched_true_cluster = true_k,
      path = names(estimated_values),
      estimated = as.numeric(estimated_values),
      truth = as.numeric(true_values),
      difference = as.numeric(difference),
      stringsAsFactors = FALSE
    )
  }

  loading_table <- do.call(rbind, loading_tables)
  rownames(loading_table) <- NULL

  # ---------------------------------------------------------------------------
  # 3. Factor covariance recovery
  # ---------------------------------------------------------------------------

  ngroups <- dim(cov_eta)[3]

  lower_index <- lower.tri(cov_eta[, , 1], diag = TRUE)
  diag_index <- diag(nrow(cov_eta[, , 1])) == 1
  offdiag_index <- lower.tri(cov_eta[, , 1], diag = FALSE)

  covariance_tables <- list()

  for (g in seq_len(ngroups)) {
    estimated_cov <- cov_eta[, , g]
    true_cov <- true_cov_eta[, , g]

    rownames(estimated_cov) <- latent_vars
    colnames(estimated_cov) <- latent_vars
    rownames(true_cov) <- latent_vars
    colnames(true_cov) <- latent_vars

    covariance_tables[[g]] <- data.frame(
      group = g,
      type = ifelse(diag_index[lower_index], "variance", "covariance"),
      path = paste(
        rownames(estimated_cov)[row(estimated_cov)[lower_index]],
        "~~",
        colnames(estimated_cov)[col(estimated_cov)[lower_index]]
      ),
      estimated = as.numeric(estimated_cov[lower_index]),
      truth = as.numeric(true_cov[lower_index]),
      difference = as.numeric(estimated_cov[lower_index] - true_cov[lower_index]),
      stringsAsFactors = FALSE
    )
  }

  covariance_table <- do.call(rbind, covariance_tables)
  rownames(covariance_table) <- NULL

  variance_table <- covariance_table[covariance_table$type == "variance", ]
  factor_cov_table <- covariance_table[covariance_table$type == "covariance", ]

  # ---------------------------------------------------------------------------
  # Return summaries
  # ---------------------------------------------------------------------------

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

    loading = list(
      paths = paths,
      table = loading_table,
      bias = mean(loading_table$difference),
      rmse = sqrt(mean(loading_table$difference^2)),
      max_abs_error = max(abs(loading_table$difference))
    ),

    factor_covariance = list(
      table = covariance_table,
      bias = mean(covariance_table$difference),
      rmse = sqrt(mean(covariance_table$difference^2)),
      max_abs_error = max(abs(covariance_table$difference)),

      variance_bias = mean(variance_table$difference),
      variance_rmse = sqrt(mean(variance_table$difference^2)),
      variance_max_abs_error = max(abs(variance_table$difference)),

      covariance_bias = mean(factor_cov_table$difference),
      covariance_rmse = sqrt(mean(factor_cov_table$difference^2)),
      covariance_max_abs_error = max(abs(factor_cov_table$difference))
    )
  )
}
