#' Run Step 1: measurement model clustering
#'
#' Fit the first step of the double-mixture multigroup SEM (2MixMG-SEM) method.
#' In this step, groups are clustered based on (a subset of) measurement model
#' parameters using mixture multigroup confirmatory factor analysis via `mixmgfa`.
#' The resulting group- and cluster-specific factor covariance matrices are weighted
#' by the posterior cluster membership probabilities to obtain one factor covariance matrix
#' for each group. These matrices are used as input to [MixMix_Step2()].
#'
#' @param data A data frame containing the observed variables and a grouping variable.
#' @param group A character string giving the name of the grouping variable. Defaults to `"group"`.
#' @param step1model A lavaan model syntax string specifying the measurement model.
#' @param MM.cluster.spec A character vector passed to `mixmgfa::mixmgfa()` specifying which measurement parameters are cluster-specific.
#'        Defaults to `"loadings"`. Use `c("loadings", "residuals")` to also make the unique variances cluster-specific.
#' @param MM.nclus Integer scalar indicating the number of measurement model clusters.
#' @param MM.maxiter Maximum number of iterations passed to `mixmgfa::mixmgfa()`. Increase in case of non-convergence.
#' @param MM.nruns Number of random starts passed to `mixmgfa::mixmgfa()` (to avoid local maxima).
#' @param MM.design A zero-one loading matrix passed to `mixmgfa::mixmgfa()`. Rows correspond to observed variables and columns to factors;
#'        `0` indicates a zero loading and `1` indicates a non-zero loading.
#' @param invar_loadings Optional invariant loading specification passed to `mixmgfa::mixmgfa()`. Defaults to `NULL`.
#' @param markers A zero-one matrix specifying the marker variable for each factor. Rows correspond to observed variables and columns to
#'        factors. An entry of `1` identifies the marker variable; all other entries should be `0`.
#' @param seed Integer seed used before fitting the measurement model.
#'
#' @return A named list containing:
#' \describe{
#'   \item{cov_eta}{An array of group-specific factor covariance matrices.}
#'   \item{ngroups}{Number of groups.}
#'   \item{N_gs}{Group sample sizes.}
#'   \item{S_unbiased}{Observed sample covariance matrices per group.}
#'   \item{vars}{Observed variable names.}
#'   \item{lat_var}{Latent variable names.}
#'   \item{mmgfa_output}{Selected outputs from the Step 1 `mixmgfa` fit.}
#'   \item{step1_time}{Step 1 computation time in minutes.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' step1_model <- '
#'   F1 =~ x1 + x2 + x3 + x4 + x5
#'   F2 =~ z1 + z2 + z3 + z4 + z5
#'   F3 =~ m1 + m2 + m3 + m4 + m5
#'   F4 =~ y1 + y2 + y3 + y4 + y5
#' '
#'
#' step1_out <- MixMix_Step1(
#'   data = data,
#'   group = "country",
#'   step1model = step1_model,
#'   MM.cluster.spec = c("loadings"),
#'   MM.nclus = 2,
#'   MM.design = design,
#'   invar_loadings = NULL,
#'   markers = markers
#' )
#' }


MixMix_Step1 <- function(data, group = "group", step1model, MM.cluster.spec = c("loadings"),
                         MM.nclus, MM.maxiter = 10000, MM.nruns = 50,
                         MM.design, invar_loadings = NULL, markers, seed = 100) {

  start_time_step1 <- Sys.time()

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (length(MM.nclus) != 1L) {
    stop(
      "This development version expects `MM.nclus` to be a single number of measurement model clusters.",
      call. = FALSE
    )
  }

  # center observed variables within groups
  g_name <- unique(data[[group]])
  partable <- lavaan::lavaanify(step1model, auto = TRUE)
  vars <- lavaan::lavNames(partable) #obs var
  lat_var <- lavaan::lavNames(partable, type = "lv") #latent var

  centered <- data
  group.idx <- match(data[[group]], g_name)
  group.sizes <- tabulate(group.idx)
  group.means <- rowsum.default(as.matrix(data[, vars, drop = FALSE]),
                                group = group.idx, reorder = FALSE,
                                na.rm = FALSE)/group.sizes
  centered[, vars] <- data[, vars, drop = FALSE] - group.means[group.idx, , drop = FALSE]

  N_gs <- group.sizes
  nfactors <- length(lat_var)
  ngroups <- length(N_gs)

  # mixmgfa expects the grouping variable in the first column
  centered <- centered[, c(group, vars), drop = FALSE]
  centered[[group]] <- group.idx

  # compute the sample covariance matrix per group
  S_unbiased <- lapply(unique(centered[[group]]),
                       function(x) {stats::cov(centered[centered[[group]] == x, vars, drop = FALSE])})

  set.seed(seed)
  # fit the MixMG-CFA model
  mixmgfa_args <- list(
    data = centered, N_gs = N_gs, nfactors = nfactors,
    cluster.spec = MM.cluster.spec, nsclust = MM.nclus,
    maxiter = MM.maxiter, nruns = MM.nruns, design = MM.design
  )

  if (!is.null(invar_loadings)) {
    if (!"invar_loadings" %in% names(formals(mixmgfa::mixmgfa))) {
      stop(
        "The installed version of `mixmgfa::mixmgfa()` does not support `invar_loadings`. ",
        "Use `invar_loadings = NULL`, or install a version that supports this argument.",
        call. = FALSE
      )
    }

    mixmgfa_args$invar_loadings <- invar_loadings # add the argument when supported
  }

  output1 <- do.call(mixmgfa::mixmgfa, mixmgfa_args)

  # rescale factors using marker variables
  output2 <- mixmgfa:::ScaleRotateMixmgfa(output1, N_gs = N_gs, cluster.spec = MM.cluster.spec,
                                nsclust = MM.nclus, design = MM.design, rescale = 1, markers = markers,
                                rotation = 0, targetT = 0, targetW = 0)
  mmgfa_solution <- output2$MMGFAsolutions[[paste0(MM.nclus, ".clusters")]]

  # unique variances can be either group-specific or cluster-specific
  uniquevar_key <- grep("uniquevariances$", names(mmgfa_solution), value = TRUE)
  # extract the parameters needed for Step 2
  factor_cov_list <-  mmgfa_solution$group.and.clusterspecific.factorcovariances
  cluster_memb <- mmgfa_solution$clustermemberships
  lambda_list <- mmgfa_solution$clusterspecific.loadings
  theta_list <- mmgfa_solution[[uniquevar_key]]

  # compute the group-specific fcov using weighted sum approach
  weighted_sum_mg <- array(data = NA_real_, dim = c(nfactors, nfactors, ngroups))

  for (g in seq_len(ngroups)) {
      # match group- and cluster-specific fcov matrices and probabilities by MM-cluster,
      # multiply and sum across MM-clusters
      weighted_sum_mg[, , g] <- Reduce(`+`, Map(`*`, factor_cov_list[g, ], cluster_memb[g, ]))
    }
  cov_eta <- weighted_sum_mg
  dimnames(cov_eta)[[1]] <- dimnames(cov_eta)[[2]] <- lat_var

  end_time_step1 <- Sys.time()
  step1_time <- difftime(end_time_step1, start_time_step1, units = "mins")

  return(list(cov_eta = cov_eta,
              ngroups = ngroups,
              N_gs = N_gs,
              S_unbiased   = S_unbiased,
              vars         = vars,
              lat_var      = lat_var,
              mmgfa_output = list(output1 = output1,
                                  output2 = output2,
                                  cluster_memb = cluster_memb,
                                  factor_cov = factor_cov_list,
                                  lambda = lambda_list,
                                  theta = theta_list),
              step1_time = step1_time))

}
